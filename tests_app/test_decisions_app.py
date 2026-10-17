"""DEC-1 (app): without --strict the engine writes probably-wrong genes as PAL2NAL does; the app must still mark
them clearly in its report and GUI table."""
import os
import tempfile
import unittest

from coal_app import core

UNREL_PROT = "MKVALWYFHIKLDEPQRSTGNCMKVALWYFHI"
UNREL_DNA = "ATGCCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT"
MASKED_PROT = "MTVDCLAHRFISDLVKSPSFGSVLACYVISPLDRNVCGQH"
MASKED_DNA = ("GCATGNNNNNNNNNNNNTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATCAGTNNNNNNNNNNNN"
              "NNNNNNTGTNNNNNNNNN")
GOOD_PROT, GOOD_DNA = "MKVAL", "ATGAAAGTTGCTCTG"


def _pad(p, n):
    return p + "-" * (n - len(p))


class TestDec1AppMarksSuspectGenes(unittest.TestCase):  # DEC-1
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_dec_app_")
        self.d = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def files(self, pep, nuc):
        p, n = os.path.join(self.d, "p.fa"), os.path.join(self.d, "n.fa")
        with open(p, "w") as fh:
            fh.write(pep)
        with open(n, "w") as fh:
            fh.write(nuc)
        return p, n

    def run_app(self, pep, nuc):
        p, n = self.files(pep, nuc)
        return core.run_coal(p, n, out_dir=os.path.join(self.d, "run"), probe_codes=False)

    def by_id(self, res):
        return {s["id"]: s for s in res.report["sequences"]}

    def test_parse_new_notices(self):
        sid, msg = core._parse_warning(
            "coal: warning: u: protein and DNA are probably unrelated: 30 of 32 checkable codons do not match the "
            "protein. The row was written as PAL2NAL writes it, but it is probably the wrong DNA or the wrong "
            "genetic code (--strict makes this an error)")
        self.assertEqual((sid, msg["code"]), ("u", "suspect_unrelated"))
        self.assertIn("30 of the 32", msg["text"])
        sid, msg = core._parse_warning(
            "coal: warning: g: position of the gene in the DNA is a guess: no 10-residue segment of the protein was "
            "found in the DNA and only 4 of 40 codons have no ambiguous bases, too few to tell where the gene "
            "starts in the longer DNA. The codons were taken from the start of the DNA as PAL2NAL does, so the row "
            "may be shifted (--strict makes this an error)")
        self.assertEqual((sid, msg["code"]), ("g", "suspect_position"))
        sid, msg = core._parse_warning(  # FZ-4 wording
            "coal: warning: s: position of the gene in the DNA is a guess: no 10-residue segment of the protein was "
            "found in the DNA and the gene has only 7 codons, too few to tell where it starts in the longer DNA. ...")
        self.assertEqual((sid, msg["code"]), ("s", "suspect_position"))
        self.assertIn("only 7 codons", msg["text"])

    def test_unrelated_gene_is_failed_but_alignment_written(self):
        n = len(UNREL_PROT)
        res = self.run_app(f">u\n{UNREL_PROT}\n>ok\n{_pad(GOOD_PROT, n)}\n", f">u\n{UNREL_DNA}\n>ok\n{GOOD_DNA}\n")
        self.assertEqual(res.exit_code, 0, res.stderr)
        self.assertIsNotNone(res.output_path)
        seqs = self.by_id(res)
        self.assertEqual(seqs["u"]["status"], "failed")
        self.assertEqual(seqs["u"]["messages"][0]["code"], "suspect_unrelated")
        self.assertEqual(seqs["ok"]["status"], "ok")
        self.assertEqual(res.report["summary"]["failed"], 1)
        sentence = core.sequence_sentence(seqs["u"])
        self.assertIn("wrong DNA", sentence)
        self.assertIn("--strict", sentence)
        rows = core.health_rows(res.report)
        self.assertTrue(any(r[2] == "u" and "Failed" in r[1] for r in rows), rows)
        self.assertNotEqual(core.banner_head(res.ok, res.report["summary"], res.hints),
                            "### ✔ Done: every sequence matched")

    def test_masked_gene_with_utr_is_failed(self):
        res = self.run_app(f">g\n{MASKED_PROT}\n", f">g\n{MASKED_DNA}\n")
        self.assertEqual(res.exit_code, 0, res.stderr)
        g = self.by_id(res)["g"]
        self.assertEqual(g["status"], "failed")
        self.assertEqual(g["messages"][0]["code"], "suspect_position")
        self.assertIn("shifted", core.sequence_sentence(g))

    def test_partial_id_match_marks_mispaired_genes(self):
        # by order, geneA gets geneB's (unrelated) DNA; geneB gets geneX's DNA, which does encode it
        n = len(UNREL_PROT)
        res = self.run_app(f">geneA\n{UNREL_PROT}\n>geneB\n{_pad(GOOD_PROT, n)}\n",
                           f">geneB\n{UNREL_DNA}\n>geneX\n{GOOD_DNA}\n")
        self.assertEqual(res.exit_code, 0, res.stderr)
        seqs = self.by_id(res)
        self.assertEqual(seqs["geneA"]["status"], "failed")
        self.assertTrue(any(m["code"] == "suspect_pairing" for m in seqs["geneA"]["messages"]))
        self.assertEqual(seqs["geneB"]["status"], "warning")
        self.assertTrue(any(m["code"] == "suspect_pairing" for m in seqs["geneB"]["messages"]))
        self.assertIn("geneX", core.sequence_sentence(seqs["geneB"]))
        self.assertTrue(any(h.get("code") == "partial_id_match" for h in res.hints))


class TestDec4AppFormats(unittest.TestCase):  # DEC-4
    def test_engine_offers_all_formats_and_app_writes_them(self):
        self.assertEqual(set(core.detect_features().output_formats), {"fasta", "clustal", "paml", "codon"})
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        fx = os.path.join(root, "tests", "fixtures", "clean")
        with tempfile.TemporaryDirectory(prefix="coal_dec_fmt_") as d:
            for fmt in ("clustal", "paml", "codon", "fasta"):
                with self.subTest(fmt=fmt):
                    res = core.run_coal(os.path.join(fx, "pep.fasta"), os.path.join(fx, "nuc.fasta"),
                                        output_format=fmt, out_dir=os.path.join(d, fmt), probe_codes=False)
                    self.assertTrue(res.ok, res.stderr)
                    with open(res.output_path, encoding="utf-8") as fh:
                        got = fh.read()
                    with open(os.path.join(root, "tests", "golden", "clean", f"{fmt}.stdout"), encoding="utf-8") as fh:
                        self.assertEqual(got, fh.read())
                    self.assertEqual(res.report["summary"]["ok"], 3)


if __name__ == "__main__":
    unittest.main()
