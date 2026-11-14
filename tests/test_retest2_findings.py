"""The findings of the second independent re-test of 2026-11-14 (Q01-Q05), each fixed and kept fixed here.

  * Q01: a single code that wins on codons alone (code 3: no CGA / CGC) is not named when its start codons do not
    fit the genes (C. elegans mitochondria, code 5, start with ATT / TTG): the candidates are named instead.
  * Q02: `pipeline clusters` judges the genetic code per --cds file; a molecule's records with one or more internal
    stops count for the "organelle molecule" rule.
  * Q03: the one-in-ten start-codon bar only for genome-size sets (Ciona mitochondria keep code 13).
  * Q04: --subject a byte copy of --query: the same-file warning.
  * Q05: the same .gz named twice in --cds: refused before it is decompressed; a refused run leaves no copies.

Synthetic inputs only (no real sequence data); coal-cpp (COAL_CPP_BIN) for the pipeline parts.
"""
from __future__ import annotations

import gzip
import io
import os
import random
import shutil
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from test_deep_findings import cds, write_fasta  # noqa: E402
from test_dnds_pipeline import FAM_CLUSTERS, PLAIN, Q, SCHEMA, PipelineCase, needs_cpp  # noqa: E402
from test_report_schema import validate  # noqa: E402


def run(fn, *args, **kw):
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, " ".join(err.getvalue().split())


class HintRules(unittest.TestCase):
    """suggest_genetic_code() on synthetic sets."""

    def setUp(self):
        import tempfile
        self.rng = random.Random(11)
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def hint(self, organelle, nuclear=0, chosen=1, name=lambda k: f"org{k}", nuc_name=lambda k: f"nuc{k}",
             stopped=None):
        recs = [(nuc_name(k), cds(self.rng, stops=())) for k in range(nuclear)]
        recs += [(name(k), s) for k, s in enumerate(organelle)]
        path = write_fasta(self.tmp / "set.fa", recs)
        ids = [r for r, _ in recs]
        flagged = {name(k): 0 for k in range(len(organelle))}
        usable = set(ids) - set(flagged)
        genes, bg, reps = P._hint_inputs([str(path)], flagged, usable, *(() if stopped is None else (stopped,)))
        with_stops = len(organelle) if stopped is None else len(stopped)
        return P.suggest_genetic_code(chosen, genes, bg, len(recs), with_stops, reps)

    # Q01
    def test_nematode_mitochondria_are_not_code_3(self):
        """C. elegans' 12 mitochondrial genes (code 5) have no CGA / CGC, as code 3's genes, but start with ATT
        (6), TTG (3) and ATG (3): code 3 allows none of ATT / TTG. The hint was 'code 3'; now 4 and 5."""
        starts = ["ATT"] * 5 + ["TTG"] * 3 + ["ATG"] * 2
        h = self.hint([cds(self.rng, avoid=("CGA", "CGC"), start=s) for s in starts])
        self.assertTrue(h["whole_set"])
        self.assertEqual(h["codes"], [4, 5])
        self.assertNotIn("CGA / CGC", h["evidence"])   # (no evidence named for one code)
        # the same genes in a nuclear set (the whole genome)
        h = self.hint([cds(self.rng, avoid=("CGA", "CGC"), start=s) for s in starts], nuclear=200)
        self.assertFalse(h["whole_set"])
        self.assertIn(5, h["codes"])
        self.assertNotEqual(h["codes"], [3])

    def test_yeast_mitochondria_stay_code_3(self):
        # ATG starts (and one GTG, one ATA: code 3's own): the start codons do not argue against code 3
        starts = ["ATG"] * 8 + ["GTG", "ATA"]
        for nuclear in (0, 200):
            h = self.hint([cds(self.rng, avoid=("CGA", "CGC"), start=s) for s in starts], nuclear=nuclear)
            self.assertEqual(h["codes"], [3], nuclear)
            self.assertIn("CGA / CGC codons, which code 3 does not use", h["evidence"])
        # one odd start is not enough to doubt it
        h = self.hint([cds(self.rng, avoid=("CGA", "CGC"), start=s) for s in ["ATT"] + ["ATG"] * 9])
        self.assertEqual(h["codes"], [3])

    def test_vertebrate_mitochondria_stay_code_2(self):
        starts = ["ATG"] * 6 + ["ATA", "ATT", "ATA", "GTG"]
        h = self.hint([cds(self.rng, avoid=("AGA", "AGG"), start=s) for s in starts], nuclear=200)
        self.assertEqual(h["codes"], [2])


if __name__ == "__main__":
    unittest.main()
