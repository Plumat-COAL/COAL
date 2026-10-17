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
        self.assertTrue(r["completed"])
        self.assertTrue(r["all_ok"])  # IMP-M-1
        self.assertFalse(r["needs_attention"])
        self.assertEqual(r["counts"], {"total": 8, "ok": 8, "warning": 0, "failed": 0, "not_checked": 0})
        self.assertIsNone(r["suggested_codontable"])
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
        # CHANGED(IMP-M-1): was {"ok": false, "error": ...}; now a tool error, so MCP clients see isError
        with self.assertRaisesRegex(tools.CoalToolError, "not found.*absolute path"):
            tools.codon_align("/no/such/file.fa", "/no/such/dna.fa")

    def test_explain_report(self):
        text = tools.explain_report(core.SAMPLE_REPORT)
        self.assertIn("Homo_COX1", text)  # APP-15
        self.assertIn("genetic code 2", text)
        with self.assertRaisesRegex(tools.CoalToolError, "Could not read"):  # CHANGED(IMP-M-1): was returned
            tools.explain_report("/no/such/report.json")

    def test_explain_after_failed_run(self):
        p, d = core.example_paths("one_broken_gene")
        r = tools.codon_align(p, d)
        self.assertIn("OG0001_sp14", tools.explain_report(r["report_path"]))

    def test_docstrings_for_agents(self):
        for fn in (tools.codon_align, tools.explain_report):
            self.assertIn("Args:", fn.__doc__)
            self.assertIn("Returns:", fn.__doc__)


class TestTruthfulResult(_PrivateRuns):
    """IMP-M-1: completed (COAL finished) is not the same as all_ok (every sequence matched)."""

    def test_wrong_code_completes_but_needs_attention(self):
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1)
        self.assertTrue(r["completed"])
        self.assertFalse(r["all_ok"])
        self.assertTrue(r["needs_attention"])
        self.assertEqual((r["counts"]["ok"], r["counts"]["warning"]), (0, 24))
        self.assertEqual(r["suggested_codontable"], 2)
        self.assertIn("codontable=2", r["next_step"])
        self.assertNotIn("ok", r)  # the old, misleading flag is gone
        self.assertNotIn("stderr", r)

    def test_result_is_compact(self):
        import json
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1)
        self.assertLessEqual(len(r["problems"]), tools.MAX_PROBLEMS)
        self.assertEqual(r["problems_total"], 24)
        self.assertGreater(r["warning_lines"], 24)
        self.assertFalse([ln for ln in r["engine_messages"] if ln.startswith("coal: warning:")])
        self.assertLess(len(json.dumps(r)), 12000)

    def test_bad_code_and_format_are_tool_errors(self):
        p, d = core.example_paths("clean")
        with self.assertRaisesRegex(tools.CoalToolError, "unknown genetic code 99"):
            tools.codon_align(p, d, codontable=99)
        with self.assertRaisesRegex(tools.CoalToolError, "output format 'xml'"):
            tools.codon_align(p, d, output_format="xml")

    def test_summarise_counts_and_hint(self):
        rep = core.load_report(core.SAMPLE_REPORT)
        res = core.RunResult(exit_code=0, output_path="/x/out.fasta", report_path=core.SAMPLE_REPORT, report=rep,
                             stderr="coal: warning: a pepAlnPos 1: M does not correspond to ATA\n"
                                    "coal: warning: a pepAlnPos 2: W does not correspond to TGA\nsome notice\n",
                             coal_cmd="coal a b", pal2nal_cmd="pal2nal.pl a b",
                             hints=[{"code": "genetic_code", "suggested_codontable": 2, "text": "t"}])
        out = tools.summarise(res)
        self.assertTrue(out["completed"])
        self.assertFalse(out["all_ok"])
        self.assertEqual(out["suggested_codontable"], 2)
        self.assertEqual((out["warning_lines"], out["engine_messages"]), (2, ["some notice"]))
        self.assertEqual(set(out), set(tools.CodonAlignResult.__annotations__))

    def test_a_genetic_code_hint_means_not_all_ok_even_if_counts_are_clean(self):
        rep = {"summary": {"total": 2, "ok": 2, "warning": 0, "failed": 0}, "sequences": []}
        res = core.RunResult(0, "/x/o.fasta", "/x/r.json", rep, "", "coal", "pal2nal.pl",
                             [{"code": "partial_id_match", "text": "t"}])
        self.assertFalse(tools.summarise(res)["all_ok"])
        res.hints = []
        self.assertTrue(tools.summarise(res)["all_ok"])


class TestToolSchemas(unittest.TestCase):
    """IMP-M-1: argument descriptions and enums (built from the type hints by mcp and by Gradio)."""

    def test_enums_match_what_coal_supports(self):
        import typing
        self.assertEqual(list(typing.get_args(tools.GeneticCode)), list(core.GENETIC_CODES))
        self.assertEqual(set(typing.get_args(tools.OutputFormat)), set(core.OUTPUT_FORMATS))

    def test_paths_are_described_as_absolute(self):
        import typing
        for fn, names in ((tools.codon_align, ("protein_alignment_path", "dna_fasta_path")),
                          (tools.explain_report, ("report_path",))):
            hints = typing.get_type_hints(fn, include_extras=True)
            for n in names:
                meta = [m for m in getattr(hints[n], "__metadata__", ()) if m is not None]
                if not meta:
                    self.skipTest("pydantic not installed: no Field metadata")
                self.assertIn("Absolute path", meta[0].description, n)

    @unittest.skipUnless(importlib.util.find_spec("pydantic"), "pydantic not installed")
    def test_json_schema_has_enums_and_descriptions(self):
        import inspect
        from pydantic import TypeAdapter
        sig = inspect.signature(tools.codon_align, eval_str=True)
        code = TypeAdapter(sig.parameters["codontable"].annotation).json_schema()
        fmt = TypeAdapter(sig.parameters["output_format"].annotation).json_schema()
        self.assertEqual(code["enum"], list(core.GENETIC_CODES))
        self.assertEqual(fmt["enum"], ["fasta", "clustal", "paml", "codon"])
        self.assertIn("vertebrate mitochondrial", code["description"])
        out = TypeAdapter(tools.CodonAlignResult).json_schema()
        for k in ("completed", "all_ok", "needs_attention", "counts", "suggested_codontable"):
            self.assertIn(k, out["required"])

    @unittest.skipUnless(importlib.util.find_spec("gradio"), "gradio not installed")
    def test_gradio_api_schema(self):
        from gradio.utils import get_function_params, python_type_to_json_schema
        params = {name: python_type_to_json_schema(t) for name, _, _, t in get_function_params(tools.codon_align)}
        self.assertEqual(params["codontable"].get("enum"), list(core.GENETIC_CODES))
        self.assertEqual(params["protein_alignment_path"].get("type"), "string")


@unittest.skipUnless(has_mcp2() or REQUIRE_MCP, "mcp>=2 not installed (pip install -e .[mcp])")
class TestStdioMcpServer(_PrivateRuns):
    def test_lists_and_calls_tools(self):
        from coal_app import mcp_stdio
        s = mcp_stdio.build_server()
        names = sorted(t.name for t in asyncio.run(s.list_tools()))
        self.assertEqual(names, ["codon_align", "explain_report"])
        res = asyncio.run(s.call_tool("explain_report", {"report_path": core.SAMPLE_REPORT}))
        self.assertIn("Homo_COX1", res.content[0].text)  # APP-15

    def test_schemas_errors_and_version(self):  # IMP-M-1
        from coal import __version__
        from coal_app import mcp_stdio
        from mcp.server.mcpserver.exceptions import ToolError
        s = mcp_stdio.build_server()
        self.assertEqual(s.version, __version__)
        t = {t.name: t for t in asyncio.run(s.list_tools())}["codon_align"]
        props = t.input_schema["properties"]
        self.assertEqual(props["codontable"]["enum"], list(core.GENETIC_CODES))
        self.assertEqual(props["output_format"]["enum"], ["fasta", "clustal", "paml", "codon"])
        self.assertIn("Absolute path", props["dna_fasta_path"]["description"])
        self.assertIn("all_ok", t.output_schema["properties"])
        p, d = core.example_paths("mito")
        for bad in ({"protein_alignment_path": "/no/such.fa", "dna_fasta_path": d},
                    {"protein_alignment_path": p, "dna_fasta_path": d, "codontable": 99},
                    {"protein_alignment_path": p, "dna_fasta_path": d, "output_format": "xml"}):
            with self.assertRaises(ToolError):
                asyncio.run(s.call_tool("codon_align", bad))
        res = asyncio.run(s.call_tool("codon_align", {"protein_alignment_path": p, "dna_fasta_path": d}))
        self.assertEqual(res.structured_content["suggested_codontable"], 2)
        self.assertFalse(res.structured_content["all_ok"])



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
            {"jsonrpc": "2.0", "id": 4, "method": "tools/call",  # IMP-M-1: a real failure is isError
             "params": {"name": "explain_report", "arguments": {"report_path": "/no/such/report.json"}}},
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
            while 4 not in replies:
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
        self.assertEqual(sorted(k for k in replies if k is not None), [1, 2, 3, 4])
        self.assertEqual(replies[1]["result"]["serverInfo"]["name"], "coal")
        from coal import __version__
        self.assertEqual(replies[1]["result"]["serverInfo"]["version"], __version__)  # IMP-M-1: was empty
        self.assertTrue(replies[4]["result"]["isError"])
        self.assertIn("Could not read the report", replies[4]["result"]["content"][0]["text"])
        self.assertFalse(replies[3]["result"].get("isError"))
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
