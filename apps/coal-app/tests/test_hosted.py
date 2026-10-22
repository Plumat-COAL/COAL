"""The hosted Space (COAL_HOSTED=1): tools are limited to a sandbox, the intro says files go to a shared server, and MCP is off there."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core  # noqa: E402
from helpers import APP, RunsDirCase, _tmp, _write  # noqa: E402


class TestHostedSandbox(RunsDirCase):
    def _hosted(self, on):
        old = os.environ.get("COAL_HOSTED")
        if on:
            os.environ["COAL_HOSTED"] = "1"
        else:
            os.environ.pop("COAL_HOSTED", None)
        self.addCleanup(lambda: os.environ.__setitem__("COAL_HOSTED", old) if old is not None
                        else os.environ.pop("COAL_HOSTED", None))

    def test_hosted_tools_refuse_paths_outside_the_sandbox(self):
        from coal_app import tools
        self._hosted(True)
        t = _tmp(self)
        secret = _write(os.path.join(t, "secret.fasta"), ">s1\nMAKL\n")
        dna = _write(os.path.join(t, "dna.fasta"), ">s1\nATGGCTAAACTG\n")
        # a tool error (isError over MCP)
        with self.assertRaisesRegex(tools.CoalToolError, "not allowed"):
            tools.codon_align(secret, dna)
        for bad in (secret, "/etc/passwd", os.path.join(core.EXAMPLES_DIR, "..", "pyproject.toml")):
            with self.assertRaisesRegex(tools.CoalToolError, "not allowed"):
                tools.explain_report(bad)
        # a symlink inside the sandbox that points outside is refused too
        link = os.path.join(core.RUNS_DIR, "link.fasta")
        try:
            os.symlink(secret, link)
        except (OSError, NotImplementedError):
            self.skipTest("symlinks not available")
        with self.assertRaisesRegex(tools.CoalToolError, "not allowed"):
            tools.codon_align(link, dna)

    def test_hosted_tools_still_work_on_examples_and_their_own_reports(self):
        from coal_app import tools
        self._hosted(True)
        p, d = core.example_paths("clean")
        r = tools.codon_align(p, d)
        self.assertTrue(r["completed"], r)
        self.assertIn("Nothing went wrong", tools.explain_report(r["report_path"]))

    def test_local_tools_accept_any_path(self):
        from coal_app import tools
        self._hosted(False)
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGCTAAACTG\n")
        self.assertTrue(tools.codon_align(p, d)["completed"])

    def test_intro_depends_on_where_it_runs(self):
        self._hosted(False)
        # the local app shows no note
        self.assertEqual(core.intro_text(), "")
        self._hosted(True)
        text = core.intro_text()
        self.assertNotIn("Nothing leaves your computer", text)
        self.assertIn("Hugging Face", text)

    def test_space_app_marks_itself_hosted_and_turns_mcp_off(self):
        with open(os.path.join(APP, "src", "coal_app", "space_app.py"), encoding="utf-8") as fh:
            src = fh.read()
        self.assertIn('os.environ["COAL_HOSTED"] = "1"', src)
        self.assertIn("mcp=False", src)
        with open(os.path.join(APP, "space", "README.md"), encoding="utf-8") as fh:
            readme = fh.read()
        self.assertNotIn("Space as an MCP server", readme)


if __name__ == "__main__":
    unittest.main()
