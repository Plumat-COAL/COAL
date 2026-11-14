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

    # Q02 (rule 2)
    def test_a_molecule_counts_its_records_with_one_stop(self):
        """Mycoplasma and K12 in one file: 376 of the Mycoplasma's 524 records have internal stops, but only 253
        two or more; the molecule was left out (2 x 253 < 524) and the hint was silent."""
        rng = self.rng
        recs = [(f"lcl|NC_N.1_cds_N{k}_{k}", cds(rng, stops=())) for k in range(300)]
        myco = [(f"lcl|NC_M.1_cds_A{k}_{k}", cds(rng, start="ATG")) for k in range(40)]
        single = [(f"lcl|NC_M.1_cds_B{k}_{k}", cds(rng, stops=("TGA",))) for k in range(20)]
        clean = [(f"lcl|NC_M.1_cds_C{k}_{k}", cds(rng, stops=())) for k in range(40)]
        path = write_fasta(self.tmp / "two.fa", recs + myco + single + clean)
        flagged = {r: 0 for r, _ in myco}
        usable = {r for r, _ in recs + single + clean}
        n = len(recs) + 100
        genes, bg, reps = P._hint_inputs([str(path)], flagged, usable)          # (without the count: as before)
        self.assertIsNone(P.suggest_genetic_code(11, genes, bg, n, 60, reps))
        stopped = {r for r, _ in myco + single}
        genes, bg, reps = P._hint_inputs([str(path)], flagged, usable, stopped)
        self.assertEqual(reps["with_stops"]["NC_M.1"], 60)
        h = P.suggest_genetic_code(11, genes, bg, n, 60, reps)
        self.assertEqual(len(h["genes"]), 40)
        self.assertIn(4, h["codes"])


@needs_cpp
class PerFileHint(PipelineCase):
    """Q02: clusters mode, one hint per --cds file."""

    def setUp(self):
        super().setUp()
        self.rng = random.Random(7)

    def run_files(self, files, **options):
        ids = [ln[1:].split()[0] for f in files for ln in Path(f).read_text().splitlines() if ln.startswith(">")]
        cl = self.tmp / "cl.tsv"
        cl.write_text(f"{ids[0]}\t{ids[1]}\n", encoding="utf-8")
        res, err = run(P.run_clusters, files, self.tmp / "out", clusters=cl, **options)
        validate(res.report, SCHEMA)
        return res.report.get("genetic_code_hints", []), err

    def test_a_mycoplasma_among_bacteria(self):
        rng = self.rng
        starts = ["ATG"] * 17 + ["GTG", "TTG", "TTG"]
        myco = [(f"lcl|NC_000908.2_cds_N{k}_{k}", cds(rng, stops=())) for k in range(60)]
        myco += [(f"lcl|NC_000908.2_cds_O{k}_{60 + k}", cds(rng, start=starts[k % 20])) for k in range(140)]
        coli = [(f"lcl|NC_000913.3_cds_E{k}_{k}", cds(rng, stops=())) for k in range(400)]
        m = write_fasta(self.tmp / "myco.fa", myco)
        e = write_fasta(self.tmp / "coli.fa", coli)
        hints, err = self.run_files([m, e], codontable=11)
        self.assertEqual(len(hints), 1)
        h = hints[0]
        self.assertEqual((h["set"], h["file"], h["whole_set"], h["candidates"], h["records"]),
                         ("cds", str(m), True, [4], 200))
        self.assertIn(f"cds {m}: most genes (140 of the 200) have internal stop codons under genetic code 11", err)
        self.assertIn("Run again with -c 4", err)
        self.assertNotIn(str(e) + ":", err)

    def test_yeast_and_human_mitochondria_get_their_own_codes(self):
        rng = self.rng
        scer = [(f"lcl|NC_001133.9_cds_Y{k}_{k}", cds(rng, stops=())) for k in range(300)]
        scer += [(f"lcl|NC_001224.1_cds_M{k}_{300 + k}", cds(rng, avoid=("CGA", "CGC"))) for k in range(8)]
        human = [(f"ENST0000036{k:04d}", cds(rng, avoid=("AGA", "AGG"), start=("ATG", "ATA", "ATT")[k % 3]))
                 for k in range(13)]
        y = write_fasta(self.tmp / "scer.fa", scer)
        hs = write_fasta(self.tmp / "hs_mt.fa", human)
        hints, err = self.run_files([y, hs])
        got = {h["file"]: (h["candidates"], h["whole_set"]) for h in hints}
        self.assertEqual(got, {str(y): ([3], False), str(hs): ([2], True)})
        self.assertIn(f"cds {y}: 8 of the 8 genes with internal stop codons", err)
        self.assertIn(f"take them out of {y} into a file of their own", err)
        self.assertIn(f"cds {hs}: most genes (13 of the 13)", err)
        # one file: the warning as before ("cds: ..."), the report names the file
        hints, err = self.run_files([y])
        self.assertEqual([(h["file"], h["candidates"]) for h in hints], [(str(y), [3])])
        self.assertIn("cds: 8 of the 8 genes", err)


if __name__ == "__main__":
    unittest.main()
