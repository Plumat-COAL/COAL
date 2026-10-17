"""FZ-1 .. FZ-6: findings of the independent differential fuzz against the real PAL2NAL v14
(overnight/integration/CHANGES.md). Expected outputs marked "PAL2NAL v14:" were produced by running the real PAL2NAL.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
TIMEOUT = 60


def run(module, args, cwd):
    env = dict(os.environ, PYTHONPATH=str(SRC))
    return subprocess.run([sys.executable, "-m", module, *args], cwd=cwd, env=env, capture_output=True,
                          timeout=TIMEOUT, stdin=subprocess.DEVNULL)


class Tmp(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_fz_")
        self.d = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, data):
        (self.d / name).write_bytes(data if isinstance(data, bytes) else data.encode())

    def coal(self, *args, module="coal"):
        r = run(module, list(args), self.d)
        return r.returncode, r.stdout.decode("latin-1"), r.stderr.decode("utf-8", "replace")


NUC = ">a\nATGAAAGTT\n>b\nATGCGT\n"
OUT = ">a\nATGAAAGTT\n>b\nATGCGT---\n"  # PAL2NAL v14 for every input below


class TestFZ1ClustalHashAndBlankLines(Tmp):  # FZ-1
    def test_hash_line_after_block_is_a_mark_line_not_a_sequence(self):
        self.write("p.aln", "CLUSTAL W\n\na    MKV\nb    MR-\n#    ###\n")
        self.write("n.fa", NUC)
        self.assertEqual(self.coal("p.aln", "n.fa")[:2], (0, OUT))
        # PAL2NAL v14: the '#' line is the mark line under the block
        rc, out, err = self.coal("p.aln", "n.fa", "-output", "clustal")
        self.assertEqual((rc, out), (0, "CLUSTAL W multiple sequence alignment\n\na             ATGAAAGTT\n"
                                        "b             ATGCGT---\n              #########\n\n"), err)

    def test_hash_line_inside_a_block_is_skipped(self):
        self.write("p.aln", "CLUSTAL W\n\na    MKV\n# c  ###\nb    MR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.aln", "n.fa", "-output", "clustal")
        self.assertEqual((rc, out), (0, "CLUSTAL W multiple sequence alignment\n\na             ATGAAAGTT\n"
                                        "b             ATGCGT---\n              #########\n\n"), err)

    def test_clustal_after_leading_blank_lines(self):
        self.write("p.aln", "\n\nCLUSTAL W\n\na    MKV\nb    MR-\n")
        self.write("n.fa", NUC)
        self.assertEqual(self.coal("p.aln", "n.fa")[:2], (0, OUT))

    def test_app_reads_the_same_ids(self):
        self.write("p.aln", "\nCLUSTAL W\n\na    MKV\nb    MR-\n#    ###\n")
        sys.path.insert(0, str(SRC))
        from coal_app.seqformat import read_ids
        self.assertEqual(read_ids(str(self.d / "p.aln")), ["a", "b"])


class TestFZ2NonAsciiHeaderText(Tmp):  # FZ-2
    def test_non_ascii_description_in_protein_header(self):
        self.write("p.fa", b">a Kr\xc3\xbcppel-like factor 4\nMKV\n>b\nMR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual((rc, out), (0, OUT), err)

    def test_non_ascii_id_is_a_clear_error(self):
        self.write("p.fa", b">\xc3\xbc1 x\nMKV\n>b\nMR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual((rc, out), (1, ""))
        self.assertIn("p.fa:1: non-ASCII character in the sequence ID", err)

    def test_non_ascii_in_sequence_still_an_error(self):
        self.write("p.fa", b">a\nMK\xc3\xbcV\n>b\nMR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual(rc, 1)
        self.assertIn("p.fa:2: non-ASCII character", err)


class TestFZ5NogapUsesTheSelectedCode(Tmp):  # FZ-5
    def test_every_table_upper_and_lower_case(self):
        sys.path.insert(0, str(SRC))
        from coal.tables import CODONS, SUPPORTED, get_code
        for t in SUPPORTED:
            code = get_code(t).code
            prot = "M" + "".join(code[c] for c in CODONS)
            want = "ATG" + "".join(c for c in CODONS if code[c] != "*")
            for case in (str.upper, str.lower):
                with self.subTest(table=t, case=case.__name__):
                    self.write("p.fa", f">s\n{prot}\n")
                    self.write("n.fa", f">s\n{case('ATG' + ''.join(CODONS))}\n")
                    rc, out, err = self.coal("p.fa", "n.fa", "-codontable", str(t), "-nogap", "-q")
                    self.assertEqual(rc, 0, err)
                    self.assertEqual(out.replace("\n", "").replace(">s", ""), case(want))

    def test_table2_examples(self):
        # table 2: TGA = Trp (kept), AGA = stop (removed). PAL2NAL v14 does the opposite (standard-code stops).
        self.write("p.fa", ">s\nMW*\n")
        self.write("n.fa", ">s\nATGTGAAGA\n")
        self.assertEqual(self.coal("p.fa", "n.fa", "-codontable", "2", "-nogap")[:2], (0, ">s\nATGTGA\n"))

    def test_rna_stop_removed(self):
        self.write("p.fa", ">s\nMK*\n")
        self.write("n.fa", ">s\naugaaauaa\n")
        self.assertEqual(self.coal("p.fa", "n.fa", "-nogap")[:2], (0, ">s\naugaaa\n"))


class TestFZ3TfastxNotation(Tmp):  # FZ-3; expected outputs: PAL2NAL v14
    CASES = {
        ("MK\\VL", "ATGAAAGGTTCTG"): (">a\nATGAAAG--GTTCTG\n>b\nATGAAAGCTGTTCTG\n",
                                      "              M   K   1   V   L\na             ATG AAA G-- GTT CTG\n"
                                      "              M   K   A   V   L\nb             ATG AAA GCT GTT CTG\n\n"),
        ("MK/VL", "ATGAAAGGCTG"): (">a\nATGAAA---GG-CTG\n>b\nATGAAAGCTGTTCTG\n",
                                   "              M   K   -   2   L\na             ATG AAA --- GG- CTG\n"
                                   "              M   K   A   V   L\nb             ATG AAA GCT GTT CTG\n\n"),
    }

    def test_backslash_and_slash(self):
        for (prot, dna), (fasta, codon) in self.CASES.items():
            self.write("p.fa", f">a\n{prot}\n>b\nMKAVL\n")
            self.write("n.fa", f">a\n{dna}\n>b\nATGAAAGCTGTTCTG\n")
            for fmt, want in (("fasta", fasta), ("codon", codon)):
                with self.subTest(prot=prot, fmt=fmt):
                    rc, out, err = self.coal("p.fa", "n.fa", "-output", fmt)
                    self.assertEqual((rc, out), (0, want), err)

    def test_slash_without_a_following_residue_is_still_an_error(self):
        self.write("p.fa", ">a\nMKVL/\n>b\nMKVLA\n")
        self.write("n.fa", ">a\nATGAAAGTTCTG\n>b\nATGAAAGTTCTGGCT\n")
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual((rc, out), (1, ""))
        self.assertIn("invalid protein alignment characters ['/']", err)


class TestFZ6AmbiguousShortPlacement(Tmp):  # FZ-6
    def test_notice_when_a_short_protein_fits_in_several_places(self):
        self.write("p.fa", ">s\nMK\n")
        self.write("n.fa", ">s\nATGAAACCATGAAA\n")  # MK fits at 0 and at 8
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual((rc, out), (0, ">s\nATGAAA\n"), err)  # PAL2NAL v14: the leftmost, ATGAAA
        self.assertIn("s: placement is ambiguous: 2 possible positions in the DNA for this 2-residue protein; "
                      "used the leftmost", err)
        self.assertEqual(self.coal("p.fa", "n.fa", "-q")[1:], (">s\nATGAAA\n", ""))

    def test_no_notice_when_unique_or_longer(self):
        self.write("p.fa", ">s\nMK\n")
        self.write("n.fa", ">s\nCCATGAAACC\n")
        self.assertEqual(self.coal("p.fa", "n.fa")[2], "")
        self.write("p.fa", ">s\nMKVA\n")
        self.write("n.fa", ">s\nATGAAAGTTGCTATGAAAGTTGCT\n")
        self.assertEqual(self.coal("p.fa", "n.fa")[2], "")

    def test_app_counts_it_as_a_warning(self):
        sys.path.insert(0, str(SRC))
        from coal_app.core import _parse_warning
        sid, msg = _parse_warning("coal: warning: s: placement is ambiguous: 2 possible positions in the DNA for "
                                  "this 2-residue protein; used the leftmost (as PAL2NAL)")
        self.assertEqual((sid, msg["code"]), ("s", "warning"))


class TestFZ4ShortGeneWording(Tmp):  # FZ-4
    def setUp(self):
        super().setUp()
        self.write("p.fa", ">s\nMKV---LAWE\n")  # 7 residues, 10 columns: one anchor, which fails (1 mismatch)
        self.write("n.fa", ">s\nATGAAACCCCTGGCTTGGGAATAAGG\n")

    def test_strict_error_does_not_say_7_of_7(self):
        rc, out, err = self.coal("p.fa", "n.fa", "--strict")
        self.assertEqual(rc, 1)
        self.assertNotIn("7 of 7", err)
        self.assertIn("the gene has only 7 codons: too few to tell where it starts in the longer DNA", err)

    def test_notice_wording(self):
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual(rc, 0, err)
        self.assertNotIn("7 of 7", err)
        self.assertIn("the gene has only 7 codons", err)


if __name__ == "__main__":
    unittest.main()
