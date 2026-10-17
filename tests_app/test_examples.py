import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "examples")
NAMES = ["clean", "mito", "messy_ids", "one_broken_gene"]


def coal(*args):
    return subprocess.run([sys.executable, "-m", "coal", *args], capture_output=True, text=True, timeout=120)


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


if __name__ == "__main__":
    unittest.main()
