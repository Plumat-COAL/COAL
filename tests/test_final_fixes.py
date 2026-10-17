"""FIN-*: last review round, engine side (overnight/integration/CHANGES.md)."""
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


class Tmp(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_fin_")
        self.d = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def coal(self, pep, nuc, *args):
        (self.d / "p.fa").write_text(pep)
        (self.d / "n.fa").write_text(nuc)
        env = dict(os.environ, PYTHONPATH=str(SRC))
        r = subprocess.run([sys.executable, "-m", "coal", *args, "p.fa", "n.fa"], cwd=self.d, env=env,
                           capture_output=True, timeout=TIMEOUT, stdin=subprocess.DEVNULL, text=True)
        return r.returncode, r.stdout, r.stderr


class TestFinC14NonStandardResidueWording(Tmp):  # FIN-C14
    PEP = ">a\nMKUVBLZJO\n"
    NUC = ">a\nATGAAATGAGTTGATCTGGAAATTTAG\n"

    def test_messages_name_the_residue_honestly(self):
        code, out, err = self.coal(self.PEP, self.NUC, "-output", "fasta")
        self.assertEqual(code, 0, err)
        self.assertEqual(out, ">a\nATGAAATGAGTTGATCTGGAAATTTAG\n")  # stdout unchanged (any codon accepted)
        self.assertNotIn("unknown AA type", err)
        for pos, what in ((3, "U (selenocysteine)"), (5, "B (D or N)"), (7, "Z (E or Q)"), (8, "J (I or L)"),
                          (9, "O (pyrrolysine)")):
            self.assertIn(f"coal: warning: a pepAlnPos {pos}: {what} is not checked; any codon is accepted\n", err)

    def test_nomismatch_still_removes_the_column(self):  # DEC-4f unchanged
        code, out, err = self.coal(">a\nMKBVL\n", ">a\nATGAAAGATGTTCTG\n", "-output", "fasta", "-nomismatch")
        self.assertEqual(out, ">a\nATGAAAGTTCTG\n")
        self.assertIn("B (D or N) is not checked", err)


MITO_PEP = ">nd2\nMNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKK\n>nd2b\nMNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKK\n"
MITO_NUC = (">nd2\nATTAATCCTCTTGCTCAGCCGGTGATCTATTCTACGATTTTCGCGGGTACGCTTATTACTGCATTATCTTCTCATTGGTTTTTCACTTGANTGGGGCTAG"
            "AAATGAATATACTCGCATTTATCCCTGTCCTGACAAAAAAG\n>nd2b\nATGAACCCCCTTGCACAGCCGGTTATCTATAGCACGATCTTTGCAGGAACTCTC"
            "ATTACCGCGCTTTCGAGCCATTGGTTTTTCACGTGAGTCGGCTTAGAAATAAACATGCTAGCATTCATCCCTGTGCTGACCAAAAAG\n")


class TestFinE1FallbackFirstMetUsesTheCode(Tmp):  # FIN-E1
    """ATT is a start codon of code 2; nd2 has one N (column 31), so it needs the anchor fallback."""

    def test_start_codon_of_the_selected_code_is_accepted(self):
        code, out, err = self.coal(MITO_PEP, MITO_NUC, "-codontable", "2", "-output", "fasta")
        self.assertEqual(code, 0, err)
        self.assertNotIn("pepAlnPos 1:", err)
        self.assertIn("nd2 pepAlnPos 31: V does not correspond to NTG", err)  # the N codon is still unchecked

    def test_nomismatch_keeps_the_first_column(self):
        code, out, err = self.coal(MITO_PEP, MITO_NUC, "-codontable", "2", "-output", "fasta", "-nomismatch")
        rows = out.split(">")[1:]
        self.assertTrue(rows[0].startswith("nd2\nATTAATCCT"), out)
        self.assertTrue(rows[1].startswith("nd2b\nATGAACCCC"), out)

    def test_same_as_the_exact_fit_path(self):
        clean = MITO_NUC.replace("TTGANTGG", "TTGAGTGG")
        _, out_clean, err_clean = self.coal(MITO_PEP, clean, "-codontable", "2", "-output", "fasta")
        self.assertNotIn("pepAlnPos 1:", err_clean)
        _, out_n, _ = self.coal(MITO_PEP, MITO_NUC, "-codontable", "2", "-output", "fasta")
        self.assertEqual(out_n.replace("ANTGG", "AGTGG"), out_clean)

    def test_not_a_start_in_code_1(self):
        _, _, err = self.coal(MITO_PEP, MITO_NUC, "-output", "fasta")
        self.assertIn("nd2 pepAlnPos 1: M does not correspond to ATT", err)


class TestFinE2BlockWriterBoundedMemory(unittest.TestCase):  # FIN-E2
    def _write(self, fmt, memory):
        import io
        import random
        import sys
        sys.path.insert(0, str(SRC))
        from coal import pipeline
        rng = random.Random(3)
        n, alilen = 37, 3 * 47
        ids = [f"s{k}" for k in range(n)]
        rows = b"".join(bytes("".join(rng.choice("ACGT-") for _ in range(alilen)), "ascii") for _ in range(n))
        aalen = alilen // 3
        aas = b"".join(bytes("".join(rng.choice("MKV-") for _ in range(aalen)), "ascii") for _ in range(n))
        mask = "".join(rng.choice("# ") for _ in range(alilen))
        out = io.StringIO()
        pipeline.write_blocks(out, fmt, ids, io.BytesIO(rows), alilen, io.BytesIO(aas) if fmt == "codon" else None,
                              aalen, mask, memory=memory)
        return out.getvalue()

    def test_chunked_path_gives_the_same_bytes(self):
        from unittest import mock
        from coal import pipeline
        for fmt in ("clustal", "codon"):
            whole = self._write(fmt, 10 ** 9)
            for memory in (1, 700, 2000):  # chunks of 1, 2-3 and 6-8 sequences
                with self.subTest(fmt=fmt, memory=memory), \
                        mock.patch.object(pipeline, "_write_blocks_chunked", wraps=pipeline._write_blocks_chunked) as w:
                    self.assertEqual(self._write(fmt, memory), whole)
                    self.assertTrue(w.called)


class TestFinE3DropInParsesLikePal2nal(unittest.TestCase):  # FIN-E3; expectations observed with the real PAL2NAL v14
    def t(self, *argv):
        import sys
        sys.path.insert(0, str(SRC))
        from coal_compat.cli import translate
        return translate(list(argv))

    def test_dangling_options_are_ignored(self):
        self.assertEqual(self.t("p", "n", "-output", "fasta", "-output"), ["-output", "fasta", "--", "p", "n"])
        self.assertEqual(self.t("p", "n", "-output", "fasta", "-codontable"), ["-output", "fasta", "--", "p", "n"])
        self.assertEqual(self.t("p", "n", "-codontable", "2", "-codontable"),
                         ["-output", "clustal", "-codontable", "2", "--", "p", "n"])

    def test_codontable_numbers_as_perl_reads_them(self):
        for v, table in (("1.0", 1), ("2x", 2), (" 2", 2), ("2e0", 2), ("+2", 2), ("2.0", 2), ("11abc", 11)):
            with self.subTest(v=v):
                self.assertEqual(self.t("p", "n", "-codontable", v),
                                 ["-output", "clustal", "-codontable", str(table), "--", "p", "n"])

    def test_invalid_numbers_are_still_errors(self):
        from coal_compat.cli import UsageError
        for v in ("2.5", "0x2", "7", "abc", "-html", "-nostderr"):
            with self.subTest(v=v), self.assertRaises(UsageError) as cm:
                self.t("p", "n", "-codontable", v)
            self.assertIn("invalid codontable number", str(cm.exception))

    def test_flags_are_read_while_a_codontable_value_is_pending(self):
        self.assertEqual(self.t("p", "n", "-output", "fasta", "-codontable", "-nogap"),
                         ["-output", "fasta", "-nogap", "--", "p", "n"])
        self.assertEqual(self.t("p", "n", "-codontable", "-nogap", "2"),
                         ["-output", "clustal", "-codontable", "2", "-nogap", "--", "p", "n"])
        from coal_compat.cli import UsageError
        with self.assertRaises(UsageError):  # PAL2NAL: "invalid codontable number, m.fa!!"
            self.t("-codontable", "-nogap", "m.fa", "mn.fa")

    def test_output_takes_the_next_argument(self):
        from coal_compat.cli import UsageError
        with self.assertRaises(UsageError):  # PAL2NAL: "valid output format: clustal, paml, fasta, or codon"
            self.t("p", "n", "-output", "-nogap")


if __name__ == "__main__":
    unittest.main()
