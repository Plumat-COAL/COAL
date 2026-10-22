"""The example inputs in examples/: present, small, LF everywhere, and the engine gives the documented result on each."""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "examples")
NAMES = ["clean", "mito", "messy_ids", "one_broken_gene"]


def coal(*args):
    return subprocess.run([sys.executable, "-m", "coal", *args], capture_output=True, text=True, timeout=120,
                          stdin=subprocess.DEVNULL)


class TestExamples(unittest.TestCase):
    def test_files_present_and_small(self):
        for n in NAMES:
            d = os.path.join(EX, n)
            for f in ("protein.aln.fasta", "dna.fasta", "README.md"):
                self.assertTrue(os.path.isfile(os.path.join(d, f)), f"{n}/{f}")
            size = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))
            self.assertLess(size, 50 * 1024, n)

    def test_clean_and_mito_code2_pass(self):
        with tempfile.TemporaryDirectory() as t:
            for n, extra in (("clean", []), ("mito", ["-codontable", "2"])):
                d = os.path.join(EX, n)
                r = coal(os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta"),
                         "-o", os.path.join(t, "out.fa"), *extra)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertNotIn("warning", r.stderr)

    def test_mito_needs_code2(self):
        d = os.path.join(EX, "mito")
        with tempfile.TemporaryDirectory() as t:
            r = coal(os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta"),
                     "-o", os.path.join(t, "out.fa"))
        self.assertIn("does not correspond", r.stderr)

    def test_broken_gene_fails(self):
        d = os.path.join(EX, "one_broken_gene")
        with tempfile.TemporaryDirectory() as t:
            r = coal(os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta"),
                     "-o", os.path.join(t, "out.fa"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("OG0001_sp14", r.stderr)


class TestExamplesAreLF(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
