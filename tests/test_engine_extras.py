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


if __name__ == "__main__":
    unittest.main()
