"""Engine extras: summary/hints, --keep-going, safe ID pairing, methods text,
--max-warnings, --progress and --resume. Inputs are built in a temp dir; every
subprocess call has a timeout."""
from __future__ import annotations

import io
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import pipeline  # noqa: E402
from test_report import run, write_genes  # noqa: E402

FORMATS = ("fasta", "clustal", "paml", "codon")


def drop_record(path: Path, idx: int):
    """Remove record number idx from a one-line-per-sequence FASTA file."""
    lines = path.read_text().splitlines(True)
    del lines[2 * idx:2 * idx + 2]
    path.write_text("".join(lines), newline="\n")


class TestSummaryHint(unittest.TestCase):
    def test_mito_genes_run_with_code_1(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 40, table=2, mito=True, seed=9)
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
            self.assertIn("coal: summary: 40 of 40 sequences have codon mismatches. They all fit genetic code 2 "
                          "(vertebrate mitochondrial): try -codontable 2.", err)
            r = json.loads((d / "r.json").read_text())
            self.assertEqual(r["hints"][0]["suggested_codontable"], 2)


class TestKeepGoing(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = Path(self._tmp.name)
        write_genes(self.d, 20, broken={7}, seed=4)
        self.ref = self.d / "ref"
        self.ref.mkdir()
        for f in ("pep.fasta", "nuc.fasta"):
            (self.ref / f).write_text((self.d / f).read_text(), newline="\n")
            drop_record(self.ref / f, 7)

    def tearDown(self):
        self._tmp.cleanup()

    def test_default_still_fails(self):
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta"])
        self.assertEqual(rc, 1)
        self.assertEqual(out, "")

    def test_one_bad_gene_in_20(self):
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--keep-going", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(rc, 2)
        self.assertEqual((self.d / "o.fa").read_text().count(">"), 19)
        self.assertIn("coal: warning: skipped g7: protein and DNA are inconsistent", err)
        self.assertIn("coal: summary: 1 sequence of 20 could not be converted and was left out (first: g7:", err)
        r = json.loads((self.d / "r.json").read_text())
        self.assertEqual(r["exit_code"], 2)
        self.assertTrue(r["options"]["keep_going"])
        self.assertEqual({k: r["summary"][k] for k in ("total", "ok", "failed")}, {"total": 20, "ok": 19, "failed": 1})
        bad = [s for s in r["sequences"] if s["status"] == "failed"]
        self.assertEqual([(s["id"], s["messages"][0]["code"]) for s in bad], [("g7", "no_match")])
        self.assertEqual(r["output"]["sequences"], 19)

    def test_output_equals_run_without_the_bad_gene(self):
        for fmt in FORMATS:
            for extra in ([], ["-nogap"], ["-nomismatch"]):
                if fmt == "codon" and extra:
                    continue
                with self.subTest(fmt=fmt, extra=extra):
                    rc, out, _ = run(self.d, ["pep.fasta", "nuc.fasta", "--keep-going", "-output", fmt, *extra])
                    self.assertEqual(rc, 2)
                    rc2, ref, _ = run(self.ref, ["pep.fasta", "nuc.fasta", "-output", fmt, *extra])
                    self.assertEqual(rc2, 0)
                    self.assertEqual(out, ref)

    def test_parallel_workers(self):
        old = pipeline.BATCH_NT
        pipeline.BATCH_NT = 500            # many small batches over 3 processes
        try:
            for fmt in FORMATS:
                outs = []
                for jobs in (1, 3):
                    buf = io.StringIO()
                    cwd = os.getcwd()
                    os.chdir(self.d)
                    try:
                        with redirect_stdout(buf), redirect_stderr(io.StringIO()):
                            n = pipeline.run("pep.fasta", ["nuc.fasta"], jobs=jobs, min_nt_per_job=0,
                                             fmt=fmt, keep_going=True)
                    finally:
                        os.chdir(cwd)
                    self.assertEqual(n, 1)
                    outs.append(buf.getvalue())
                self.assertEqual(outs[0], outs[1], fmt)
        finally:
            pipeline.BATCH_NT = old

    def test_all_failed_is_an_error(self):
        write_genes(self.d, 3, broken={0, 1, 2})
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--keep-going"])
        self.assertEqual(rc, 1)
        self.assertEqual(out, "")
        self.assertIn("none of the 3 sequences could be converted", err)



class TestIdPairing(unittest.TestCase):
    """Protein IDs g0..g4; DNA headers renamed / reordered."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = Path(self._tmp.name)
        write_genes(self.d, 5, seed=6)
        rc, self.ref, _ = run(self.d, ["pep.fasta", "nuc.fasta"])
        self.assertEqual(rc, 0)
        self.recs = (self.d / "nuc.fasta").read_text().split(">")[1:]

    def tearDown(self):
        self._tmp.cleanup()

    def write_dna(self, names, order=(0, 1, 2, 3, 4)):
        recs = [names[i] + "\n" + self.recs[i].split("\n", 1)[1] for i in range(5)]
        (self.d / "nuc.fasta").write_text("".join(">" + recs[i] for i in order), newline="\n")

    def pairing(self):
        r = json.loads((self.d / "r.json").read_text())
        return r["summary"]["pairing"], [s["pairing"] for s in r["sequences"]], [s["dna_id"] for s in r["sequences"]]

    def check_normalised(self, names, order):
        self.write_dna(names, order)
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(rc, 0, err)
        self.assertEqual((self.d / "o.fa").read_text(), self.ref)
        self.assertNotIn("paired by file order", err)
        mode, per, dna_ids = self.pairing()
        self.assertEqual(mode, "normalised")
        self.assertEqual(dna_ids, [names[i] for i in range(5)])
        self.assertIn("IDs matched only after normalising", err)
        return per

    def test_version_suffixes_shuffled(self):
        per = self.check_normalised([f"g{i}.{i + 1}" for i in range(5)], (4, 2, 0, 1, 3))
        self.assertEqual(per, ["normalised"] * 5)

    def test_pipe_headers(self):
        self.check_normalised([f"sp|g{i}|G{i}_HUMAN" for i in range(5)], (1, 0, 2, 3, 4))
        self.check_normalised([f"tr|P0000{i}|g{i}" for i in range(5)], (3, 1, 2, 0, 4))

    def test_case_insensitive_and_mixed_exact(self):
        per = self.check_normalised(["g0", "G1", "g2", "G3.1", "g4"], (2, 3, 4, 0, 1))
        self.assertEqual(per, ["id", "normalised", "id", "normalised", "id"])

    def test_ambiguous_falls_back_to_order_with_warning(self):
        # g1.1 and g1.2 both normalise to g1: not one-to-one -> file order, with a clear warning
        self.write_dna(["g0", "g1.1", "g1.2", "g3", "g4"])
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertIn("coal: warning: IDs don't match, so sequences were paired by file order "
                      "(protein g0 <-> DNA g0, g1 <-> g1.1, ...). Check this is what you want.", err)
        self.assertEqual(self.pairing()[0], "order")
        self.assertEqual((self.d / "o.fa").read_text(), self.ref)   # same order, so same result

    def test_unrelated_ids_use_order(self):
        self.write_dna([f"n{i}" for i in range(5)])
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--report", "r.json"])
        self.assertEqual(rc, 0)
        self.assertEqual(out, self.ref)
        self.assertIn("paired by file order (protein g0 <-> DNA n0, g1 <-> n1, ...)", err)
        mode, per, dna_ids = self.pairing()
        self.assertEqual((mode, per[0], dna_ids[0]), ("order", "order", "n0"))



class TestMethodsAndMaxWarnings(unittest.TestCase):
    def test_methods_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 12, table=2, seed=2)
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "-codontable", "2", "-output", "paml", "-o", "o.paml",
                                 "--report", "r.json"])
            self.assertEqual(rc, 0, err)
            from coal import __version__
            m = json.loads((d / "r.json").read_text())["methods_text"]
            self.assertEqual(m, f"Codon alignments were generated with COAL v{__version__} "
                                "(https://github.com/chloe-newt/COAL) from 12 aligned protein sequences and their "
                                "coding sequences, using the vertebrate mitochondrial genetic code (NCBI table 2) and "
                                "PAML output; outputs are verified byte-identical to PAL2NAL v14 on the COAL test suite.")

    def test_methods_text_keep_going_and_nogap(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 5, broken={1})
            run(d, ["pep.fasta", "nuc.fasta", "--keep-going", "-nogap", "-o", "o.fa", "--report", "r.json"])
            m = json.loads((d / "r.json").read_text())["methods_text"]
            self.assertIn("from 4 aligned protein sequences and their coding sequences (1 sequence that could not "
                          "be matched to its DNA was excluded), using the standard genetic code (NCBI table 1)", m)
            self.assertIn("removing codon columns containing gaps or stop codons (-nogap)", m)

    def test_max_warnings(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 6, table=2, mito=True)
            rc, _, full = run(d, ["pep.fasta", "nuc.fasta"])
            n = full.count("does not correspond")
            self.assertGreater(n, 5)
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "--max-warnings", "5"])
            lines = [l for l in err.splitlines() if l.startswith("coal: warning:")]
            self.assertEqual(lines[:5], [l for l in full.splitlines() if "does not correspond" in l][:5])
            self.assertEqual(lines[5], f"coal: warning: ... and {n - 5} more warnings "
                                       "(see them all with --report FILE.json)")
            self.assertEqual(len(lines), 6)
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "--max-warnings", "5", "-o", "o.fa", "--report", "r.json"])
            self.assertIn(f"... and {n - 5} more warnings (see the report)", err)
            r = json.loads((d / "r.json").read_text())   # the report keeps every warning
            self.assertEqual(sum(len(s["messages"]) for s in r["sequences"]), n)
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "--max-warnings", "0"])
            self.assertEqual(err.count("does not correspond"), n)    # 0 = unlimited


if __name__ == "__main__":
    unittest.main()
