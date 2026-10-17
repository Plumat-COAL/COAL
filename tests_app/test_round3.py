"""Round 3 review fixes (overnight/integration/CHANGES.md, R3-A*, R3-D*, R3-E*)."""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from coal_app import core, tools

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FX = os.path.join(ROOT, "tests", "fixtures")


class RunsDirCase(unittest.TestCase):
    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def tmp(self):
        d = tempfile.mkdtemp(prefix="coal_test_r3_")
        self.addCleanup(shutil.rmtree, d, True)
        return d


class TestR3D1ExamplesAreLF(unittest.TestCase):  # R3-D1
    def test_gitattributes_forces_lf_for_examples(self):
        with open(os.path.join(ROOT, ".gitattributes"), encoding="utf-8") as fh:
            rules = [ln.split() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
        self.assertIn(["examples/**", "text", "eol=lf"], rules)

    def test_git_reports_lf_for_example_files(self):
        if not os.path.isdir(os.path.join(ROOT, ".git")) and not os.path.isfile(os.path.join(ROOT, ".git")):
            self.skipTest("not a git checkout")
        if shutil.which("git") is None:
            self.skipTest("git not installed")
        files = sorted(glob.glob(os.path.join(ROOT, "examples", "*", "*.fasta")))
        self.assertTrue(files)
        rel = [os.path.relpath(f, ROOT).replace(os.sep, "/") for f in files]
        out = subprocess.run(["git", "check-attr", "eol", "--"] + rel, cwd=ROOT, capture_output=True,
                             text=True, timeout=60).stdout
        for line in out.splitlines():
            self.assertTrue(line.endswith(": eol: lf"), line)

    def test_example_files_have_no_carriage_returns(self):
        files = glob.glob(os.path.join(ROOT, "examples", "*", "*.fasta"))
        self.assertTrue(files)
        for f in files:
            with open(f, "rb") as fh:
                self.assertNotIn(b"\r", fh.read(), f)


class TestR3D2CiConditions(unittest.TestCase):  # R3-D2
    def setUp(self):
        import re
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
            self.ci = fh.read()
        self.re = re

    def _job(self, name):
        m = self.re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.ci, self.re.S | self.re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def _step(self, job, name):
        m = self.re.search(rf"^(\s*)- name: {self.re.escape(name)}\n(.*?)(?=^\1- name:|\Z)", self._job(job),
                           self.re.S | self.re.M)
        self.assertIsNotNone(m, name)
        return m.group(2)

    def test_desktop_step_runs_even_after_app_tests_fail(self):
        step = self._step("python", "Desktop wrapper tests")
        cond = self.re.search(r"^\s*if: (.*)$", step, self.re.M)
        self.assertIsNotNone(cond)
        self.assertIn("!cancelled()", cond.group(1))
        self.assertIn("matrix.os != 'macos-latest'", cond.group(1))
        self.assertNotIn("always()", cond.group(1))

    def test_explainer_job_is_required(self):
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._job("explainer"), self.re.M))

    def test_windows_gui_job_still_informational(self):
        self.assertIn("continue-on-error: ${{ matrix.os == 'windows-latest' }}", self._job("gui"))


class TestR3E1AppParsesNewNotice(unittest.TestCase):  # R3-E1
    def test_new_notice_with_ambiguous_count(self):
        line = ("coal: warning: g: no 10-residue segment of the protein was found in the DNA; codons were assigned "
                "by position (1 of 12 checkable codons do not match; 30 contain N or other ambiguous bases)")
        sid, msg = core._parse_warning(line)
        self.assertEqual((sid, msg["code"]), ("g", "no_anchor"))
        self.assertIn("1 of 12 checkable codons do not match; 30 contain N", msg["text"])


if __name__ == "__main__":
    unittest.main()
