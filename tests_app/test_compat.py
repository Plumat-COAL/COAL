import os
import subprocess
import sys
import unittest

from coal_compat.cli import UsageError, translate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "examples", "clean")


class TestTranslate(unittest.TestCase):
    def test_options_any_order(self):
        self.assertEqual(translate(["-nogap", "p.aln", "-codontable", "2", "n.fa", "-output", "fasta"]),
                         ["p.aln", "n.fa", "-codontable", "2", "-nogap"])

    def test_two_dna_files(self):
        self.assertEqual(translate(["p", "a.fa", "b.fa", "-output", "fasta"]), ["p", "a.fa", "b.fa"])

    def test_other_formats_when_supported(self):
        self.assertEqual(translate(["p", "n", "-output", "paml"], ("fasta", "paml")), ["p", "n", "-output", "paml"])

    def test_clear_errors(self):
        for argv, text in ((["p", "n", "-blockonly", "-output", "fasta"], "blockonly"),
                           (["p", "n", "-html"], "html"),
                           (["p", "n", "-foo"], "unknown option"),
                           (["p", "n", "-output", "xml"], "unknown -output"),
                           (["p", "n", "-codontable", "x", "-output", "fasta"], "number"),
                           (["p", "-output", "fasta"], "at least one DNA"),
                           (["p", "n"], "add '-output fasta'")):
            with self.assertRaises(UsageError) as cm:
                translate(argv, ("fasta",))
            self.assertIn(text, str(cm.exception))


class TestCommand(unittest.TestCase):
    def test_runs_coal_to_stdout(self):
        r = subprocess.run([sys.executable, "-m", "coal_compat", os.path.join(EX, "protein.aln.fasta"),
                            os.path.join(EX, "dna.fasta"), "-output", "fasta"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(r.stdout.startswith(">gene01"))

    def test_help(self):
        r = subprocess.run([sys.executable, "-m", "coal_compat", "-h"], capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0)
        self.assertIn("-codontable", r.stdout)


if __name__ == "__main__":
    unittest.main()
