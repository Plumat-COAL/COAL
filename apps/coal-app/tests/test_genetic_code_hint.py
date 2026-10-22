"""The genetic-code hint: when it is given, what it says it tested, and which code it suggests."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core  # noqa: E402
from helpers import ExampleVariantCase, _recode  # noqa: E402


def _hint(r):
    return next((h for h in r.hints if h.get("code") == "genetic_code"), None)


class TestGeneticCodeHint(ExampleVariantCase):
    def test_pairing_only_warnings_get_no_code_hint(self):
        def ren(recs):
            recs[2][0], recs[4][0] = "geneX3", "geneX5"
            return recs
        r = self.run_app(*self.variant(dna=ren), probe=True)
        self.assertIsNone(_hint(r), r.hints)

    def test_hint_counts_what_was_tested(self):  # gene08 uses AGA for Arg, a stop in code 2
        def mito(g, i, c):  # read as vertebrate mitochondrial (code 2) genes
            return {"TGG": "TGA", "AGA": "CGT", "AGG": "CGT", "ATA": "ATT"}.get(c, "ATA" if c == "ATG" and i else c)
        tga = _recode(mito, genes=range(7))
        aga = _recode(lambda g, i, c: "AGA" if c in ("CGT", "CGC", "CGA", "CGG") else c, genes={7})
        r = self.run_app(*self.variant(dna=lambda recs: aga(tga(recs))), probe=True)
        h = _hint(r)
        self.assertIsNotNone(h, r.hints)
        self.assertEqual(h["suggested_codontable"], 2)
        self.assertNotIn("All", h["text"])
        self.assertIn("7 of 8 sequences match their protein", h["text"])
        self.assertNotIn("compared with", h["text"])

    def test_mito_example_wording(self):
        p, d = core.example_paths("mito")
        r = self.run_app(p, d, probe=True)
        h = _hint(r)
        self.assertEqual(h["suggested_codontable"], 2)
        self.assertIn("Under genetic code 2 (Vertebrate mitochondrial), all 24 sequences match their protein; "
                      "under code 1, none do.", h["text"])

    def test_alternative_start_codon_prefers_the_closest_code(self):
        ata = _recode(lambda g, i, c: "ATA" if i == 0 else c, genes={1, 4})
        r = self.run_app(*self.variant(dna=ata), probe=True)
        h = _hint(r)
        self.assertIsNotNone(h, r.hints)
        self.assertEqual(h["suggested_codontable"], 11, h["text"])
        self.assertIn("start codon", h["text"])
        self.assertNotIn("Mycoplasma", h["text"].split(" Codes")[0].split(" Other")[0])

    def test_frameshift_near_the_end_is_not_a_code_problem(self):
        def fs(recs):
            s = recs[2][1]
            k = len(s) - 3 * 6  # delete one base 5 codons before the stop
            recs[2][1] = s[:k] + s[k + 1:]
            return recs
        r = self.run_app(*self.variant(dna=fs), probe=True)
        self.assertEqual(r.exit_code, 0, r.stderr)
        g = next(x for x in r.report["sequences"] if x["id"] == "gene03")
        self.assertEqual(g["status"], "warning")
        text = core.sequence_sentence(g)
        self.assertIn("frameshift", text)
        self.assertNotIn("check the genetic code", text)
        self.assertIsNone(_hint(r), r.hints)


if __name__ == "__main__":
    unittest.main()
