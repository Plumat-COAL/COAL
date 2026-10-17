import asyncio
import importlib.util
import os
import unittest

from coal_app import core, tools


def has_mcp2():
    return importlib.util.find_spec("mcp") is not None and importlib.util.find_spec("mcp.server.mcpserver") is not None


class TestTools(unittest.TestCase):
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
        self.assertIn("OG0001_sp14", text)
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


@unittest.skipUnless(has_mcp2(), "mcp>=2 not installed")
class TestStdioMcpServer(unittest.TestCase):
    def test_lists_and_calls_tools(self):
        from coal_app import mcp_stdio
        s = mcp_stdio.build_server()
        names = sorted(t.name for t in asyncio.run(s.list_tools()))
        self.assertEqual(names, ["codon_align", "explain_report"])
        res = asyncio.run(s.call_tool("explain_report", {"report_path": core.SAMPLE_REPORT}))
        self.assertIn("OG0001_sp14", res.content[0].text)


if __name__ == "__main__":
    unittest.main()
