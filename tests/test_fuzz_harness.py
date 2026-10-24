"""FIN-C10: the fuzz harness's check of deliberate differences (no PAL2NAL needed)."""
from __future__ import annotations

import os
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fuzz_vs_pal2nal as fz  # noqa: E402


def _case(opts, features):
    c = fz.Case()
    c.opts, c.features = list(opts), set(features)
    return c


class TestCheckDeliberate(unittest.TestCase):
    # Coordinator's case 440 (seed 4242), reduced: a masked gene placed one base off the generator's frame, so the
    # reference's substituted IUPAC reading ('S' -> 'g', marked in lower case) straddles two codons: NGA|SAA vs Nga|gAA.
    REF = ">a\nATGNgagAATAA\n>b\nATGCGAGAATAA\n"          # PAL2NAL on the neutralised copy, without -nogap
    BASE = ">a\nATGNGASAATAA\n>b\nATGCGAGAATAA\n"         # COAL without -nogap
    NOGAP = ">a\nATGNGASAA\n>b\nATGCGAGAA\n"              # COAL -codontable 13 -nogap (TAA removed)
    OPTS = ["-codontable", "13", "-nogap"]

    def check(self, out, base=BASE, ref=REF):
        return fz.check_deliberate(_case(self.OPTS, {"iupac", "nogap_code"}), "fasta", 0, out, ref,
                                   lambda: (0, base))

    def test_case_440_is_explained(self):
        self.assertIsNone(self.check(self.NOGAP))

    def test_wrong_letter_is_still_caught(self):
        self.assertIsNotNone(self.check(self.NOGAP, base=">a\nATGNGTSAATAA\n>b\nATGCGAGAATAA\n"))

    def test_wrong_nogap_columns_are_caught(self):
        self.assertIsNotNone(self.check(">a\nATGNGASAATAA\n>b\nATGCGAGAATAA\n"))  # stop column kept
        self.assertIsNotNone(self.check(">a\nATGSAA\n>b\nATGGAA\n"))              # a non-stop column removed

    def test_iupac_all_stop_triplet_is_removed_by_nogap(self):  # case 2608: TRA = TAA/TGA, all stops in code 1
        ref = ">a\nATGtaaAAA\n"
        base = ">a\nATGTRAAAA\n"
        c = _case(["-nogap"], {"iupac"})
        self.assertIsNone(fz.check_deliberate(c, "fasta", 0, ">a\nATGAAA\n", ref, lambda: (0, base)))
        self.assertIsNotNone(fz.check_deliberate(c, "fasta", 0, base, ref, lambda: (0, base)))

    def test_substitution_must_be_the_documented_one(self):
        c = _case(["-codontable", "10"], {"genetic_code"})
        self.assertIsNone(fz.check_deliberate(c, "fasta", 0, ">a\nATGTGA\n", ">a\nATGtgc\n"))  # table 10 Cys
        self.assertIsNotNone(fz.check_deliberate(c, "fasta", 0, ">a\nATGTGG\n", ">a\nATGtgc\n"))
        # the gene placed 1 nt off the generator's frame (fuzz seed 20261024 case 1957, seed 20261025 case 1380:
        # both programs assign a mostly-N gene by position): the neutralised TGA -> tgc then straddles two codons of
        # the alignment grid, and is still the documented substitution
        self.assertIsNone(fz.check_deliberate(c, "fasta", 0, ">a\nNNTGANNNN\n", ">a\nNNtgcNNNN\n"))
        self.assertIsNone(fz.check_deliberate(c, "fasta", 0, ">a\nNTGANNNNN\n", ">a\nNtgcNNNNN\n"))
        self.assertIsNotNone(fz.check_deliberate(c, "fasta", 0, ">a\nNNTGGNNNN\n", ">a\nNNtgcNNNN\n"))  # not a TGA
        self.assertIsNotNone(fz.check_deliberate(c, "fasta", 0, ">a\nNNTGANNNA\n", ">a\nNNtgcNNNN\n"))  # and a real difference
        c = _case(["-codontable", "2"], {"first_met_code"})
        self.assertIsNone(fz.check_deliberate(c, "fasta", 0, ">a\nATTAAA\n", ">a\natgAAA\n"))  # code-2 start
        self.assertIsNotNone(fz.check_deliberate(c, "fasta", 0, ">a\nCCCAAA\n", ">a\natgAAA\n"))

    def test_parse_output_formats(self):
        self.assertEqual(fz.parse_output("paml", "   2      6\na\nATGAAA\nb\nATG---\n"),
                         [("a", "ATGAAA"), ("b", "ATG---")])
        self.assertEqual(fz.parse_output("codon", "      M   K\na     ATG AAA\n      M   -\nbb    ATG ---\n\n"),
                         [("a", "ATGAAA"), ("bb", "ATG---")])
        self.assertEqual(fz.parse_output("clustal", "CLUSTAL W multiple sequence alignment\n\na     ATG\nbb    ATG\n"
                                                    "      ###\n\na     AAA\nbb    ---\n\n"),
                         [("a", "ATGAAA"), ("bb", "ATG---")])


class TestEngineSelection(unittest.TestCase):
    """COAL_FUZZ_ENGINE (nightly.yml runs the fuzz once per engine): python in-process by default, coal-cpp with
    the PAL2NAL command line translated by the drop-in."""

    def tearDown(self):
        fz.select_engine("python")

    def test_python_is_the_default_and_runs_in_process(self):
        self.assertEqual(fz.select_engine("python"), "python")
        self.assertIs(fz.run_coal, fz._python_run_coal)
        env = dict(os.environ)
        env.pop("COAL_FUZZ_ENGINE", None)
        with unittest.mock.patch.dict(os.environ, env, clear=True):
            self.assertEqual(fz.select_engine(), "python")

    def test_cpp_without_a_binary_is_an_error_not_a_silent_python_run(self):
        with unittest.mock.patch.dict(os.environ, {"COAL_CPP_BIN": str(Path(__file__).parent / "no-such-binary")}):
            with self.assertRaises(FileNotFoundError):
                fz.select_engine("cpp")
        with self.assertRaises(ValueError):
            fz.select_engine("rust")

    def test_cpp_runs_the_binary_with_the_translated_command_line(self):
        binary = fz.cpp_binary()
        if binary is None:
            self.skipTest("coal-cpp not built (build/coal-cpp or COAL_CPP_BIN)")
        self.assertEqual(fz.select_engine("cpp"), "cpp")
        self.assertIsNot(fz.run_coal, fz._python_run_coal)
        with tempfile.TemporaryDirectory() as d:
            Path(d, "pep.fasta").write_text(">a\nMK\n")
            Path(d, "nuc.fasta").write_text(">a\nATGAAA\n")
            args = ["pep.fasta", "nuc.fasta", "-output", "fasta"]
            rc, out, err = fz.run_coal(d, args)
            self.assertEqual((rc, out), (0, ">a\nATGAAA\n"), err)
            self.assertEqual(fz._python_run_coal(d, args)[:2], (0, ">a\nATGAAA\n"))
            rc, out, err = fz.run_coal(d, ["pep.fasta", "nuc.fasta", "-output", "nonsense"])
            self.assertEqual((rc, out), (2, ""))  # a usage error of the drop-in's translation, as in-process


if __name__ == "__main__":
    unittest.main()
