import asyncio
import importlib.util
import os
import shutil
import tempfile
import unittest

from coal_app import core, tools


def has_mcp2():
    return importlib.util.find_spec("mcp") is not None and importlib.util.find_spec("mcp.server.mcpserver") is not None


# APP-16: in the CI job that installs .[mcp], COAL_REQUIRE_EXTRAS=mcp makes these tests run (and fail) instead
# of skipping when the package is missing.
REQUIRE_MCP = "mcp" in os.environ.get("COAL_REQUIRE_EXTRAS", "").split(",")


class _PrivateRuns(unittest.TestCase):
    """Tool runs write into core.RUNS_DIR: use a private one per test and delete it afterwards."""

    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)


class TestTools(_PrivateRuns):
    def test_codon_align_clean(self):
        p, d = core.example_paths("clean")
        r = tools.codon_align(p, d)
        self.assertTrue(r["ok"])
        self.assertTrue(os.path.isfile(r["output_path"]))
        self.assertEqual(r["problems"], [])
        self.assertIn("pal2nal.pl", r["pal2nal_command"])

    def test_codon_align_mito_hint_then_rerun(self):
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d)
        code = next(h["suggested_codontable"] for h in r["hints"] if h["code"] == "genetic_code")
        r2 = tools.codon_align(p, d, codontable=code)
        self.assertEqual(r2["problems"], [])

    def test_codon_align_missing_file(self):
        r = tools.codon_align("/no/such/file.fa", "/no/such/dna.fa")
        self.assertFalse(r["ok"])
        self.assertIn("not found", r["error"])

    def test_explain_report(self):
        text = tools.explain_report(core.SAMPLE_REPORT)
        self.assertIn("Homo_COX1", text)  # APP-15
        self.assertIn("genetic code 2", text)
        self.assertIn("Could not read", tools.explain_report("/no/such/report.json"))

    def test_explain_after_failed_run(self):
        p, d = core.example_paths("one_broken_gene")
        r = tools.codon_align(p, d)
        self.assertIn("OG0001_sp14", tools.explain_report(r["report_path"]))

    def test_docstrings_for_agents(self):
        for fn in (tools.codon_align, tools.explain_report):
            self.assertIn("Args:", fn.__doc__)
            self.assertIn("Returns:", fn.__doc__)


@unittest.skipUnless(has_mcp2() or REQUIRE_MCP, "mcp>=2 not installed (pip install -e .[mcp])")
class TestStdioMcpServer(_PrivateRuns):
    def test_lists_and_calls_tools(self):
        from coal_app import mcp_stdio
        s = mcp_stdio.build_server()
        names = sorted(t.name for t in asyncio.run(s.list_tools()))
        self.assertEqual(names, ["codon_align", "explain_report"])
        res = asyncio.run(s.call_tool("explain_report", {"report_path": core.SAMPLE_REPORT}))
        self.assertIn("Homo_COX1", res.content[0].text)  # APP-15



@unittest.skipUnless(has_mcp2() or REQUIRE_MCP, "mcp>=2 not installed (pip install -e .[mcp])")
class TestStdioMcpHandshake(_PrivateRuns):
    """APP-16: talk to `python -m coal_app --mcp` over real stdio, as Claude Code / Desktop do."""

    def test_initialize_list_and_call_over_stdio(self):
        import json
        import subprocess
        import sys
        msgs = [
            {"jsonrpc": "2.0", "id": 1, "method": "initialize",
             "params": {"protocolVersion": "2025-06-18", "capabilities": {},
                        "clientInfo": {"name": "coal-test", "version": "0"}}},
            {"jsonrpc": "2.0", "method": "notifications/initialized"},
            {"jsonrpc": "2.0", "id": 2, "method": "tools/list"},
            {"jsonrpc": "2.0", "id": 3, "method": "tools/call",
             "params": {"name": "explain_report", "arguments": {"report_path": core.SAMPLE_REPORT}}},
        ]
        env = dict(os.environ, COAL_RUNS_DIR=core.RUNS_DIR)
        import threading
        p = subprocess.Popen([sys.executable, "-m", "coal_app", "--mcp"], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, encoding="utf-8",
                             env=env)
        watchdog = threading.Timer(60, p.kill)  # never hang the suite
        watchdog.start()
        replies = {}
        try:
            for m in msgs:
                p.stdin.write(json.dumps(m) + "\n")
            p.stdin.flush()
            while 3 not in replies:
                line = p.stdout.readline()
                if not line:
                    break
                if line.strip():
                    obj = json.loads(line)
                    replies[obj.get("id")] = obj
        finally:
            watchdog.cancel()
            p.stdin.close()
            try:
                p.wait(timeout=20)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait(timeout=20)
            p.stdout.close()
        self.assertEqual(sorted(k for k in replies if k is not None), [1, 2, 3])
        self.assertEqual(replies[1]["result"]["serverInfo"]["name"], "coal")
        self.assertEqual(sorted(t["name"] for t in replies[2]["result"]["tools"]), ["codon_align", "explain_report"])
        text = "".join(c.get("text", "") for c in replies[3]["result"]["content"])
        self.assertIn("Homo_COX1", text)


class TestMcpDemoScript(_PrivateRuns):
    def test_demo_outputs_match_current_tools(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(root, "demo", "mcp_demo.md"), encoding="utf-8") as fh:
            md = fh.read()
        p, d = core.example_paths("one_broken_gene")
        r = tools.codon_align(p, d)
        self.assertIn(r["problems"][0]["explanation"], md)
        self.assertIn("claude mcp add --transport stdio coal -- python -m coal_app --mcp", md)


if __name__ == "__main__":
    unittest.main()
