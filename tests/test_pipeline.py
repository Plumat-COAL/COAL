"""Streaming / multi-process pipeline: results must not depend on -j or batch size."""
from __future__ import annotations

import io
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks"))
from coal import pipeline  # noqa: E402
from coal.core import CoalError  # noqa: E402
from gen_data import generate  # noqa: E402

FIX = ROOT / "tests" / "fixtures"
CASES = json.loads((ROOT / "tests" / "cases.json").read_text())["cases"]


def run(case_dir, pep, nucs, **kw):
    out = io.StringIO()
    cwd = os.getcwd()
    os.chdir(case_dir)
    try:
        with redirect_stdout(out):
            pipeline.run(pep, nucs, **kw)
    finally:
        os.chdir(cwd)
    return out.getvalue()


class TestEffectiveJobs(unittest.TestCase):
    def test_small_inputs_use_one_process(self):
        self.assertEqual(pipeline.effective_jobs(0, 1_000_000), 1)
        self.assertEqual(pipeline.effective_jobs(8, 1_000_000), 1)

    def test_large_inputs_scale_up_to_the_limit(self):
        cores = os.cpu_count() or 1
        self.assertEqual(pipeline.effective_jobs(0, 10**12), cores)
        self.assertEqual(pipeline.effective_jobs(3, 10**12), 3)
        self.assertEqual(pipeline.effective_jobs(0, 3 * pipeline.MIN_NT_PER_JOB), min(3, cores))

    def test_one_means_one(self):
        self.assertEqual(pipeline.effective_jobs(1, 10**12), 1)


class TestParallelEqualsSerial(unittest.TestCase):
    CASES = ["clean", "utr_frames", "frameshift_mixed", "long_mismatch", "two_nuc_files",
             "unordered_ids", "id_mismatch_same_order", "clustal_input"]

    def test_fixtures_any_jobs_any_batch(self):
        old = pipeline.BATCH_NT
        try:
            for case in self.CASES:
                spec = CASES[case]
                args = (FIX / case, spec["pep"], spec.get("nuc", ["nuc.fasta"]))
                for opts in ({}, {"nogap": True}, {"nomismatch": True},
                             {"fmt": "clustal"}, {"fmt": "paml"}, {"fmt": "codon"},
                             {"fmt": "paml", "nogap": True}, {"fmt": "clustal", "nomismatch": True}):
                    pipeline.BATCH_NT = 2_000_000
                    ref = run(*args, jobs=1, **opts)
                    pipeline.BATCH_NT = 1           # one sequence per batch: stresses ordering
                    with self.subTest(case=case, opts=opts):
                        self.assertEqual(run(*args, jobs=3, min_nt_per_job=0, **opts), ref)
        finally:
            pipeline.BATCH_NT = old

    def test_generated_dataset(self):
        with tempfile.TemporaryDirectory() as d:
            generate(300, 200, Path(d), seed=7)
            ref = run(d, "pep.fasta", ["nuc.fasta"], jobs=1)
            old = pipeline.BATCH_NT
            pipeline.BATCH_NT = 5000
            try:
                self.assertEqual(run(d, "pep.fasta", ["nuc.fasta"], jobs=2, min_nt_per_job=0), ref)
            finally:
                pipeline.BATCH_NT = old
            self.assertEqual(ref.count(">"), 300)
            for fmt in ("clustal", "paml", "codon"):
                ref = run(d, "pep.fasta", ["nuc.fasta"], jobs=1, fmt=fmt)
                pipeline.BATCH_NT = 5000
                try:
                    with self.subTest(fmt=fmt):
                        self.assertEqual(run(d, "pep.fasta", ["nuc.fasta"], jobs=2, min_nt_per_job=0, fmt=fmt), ref)
                finally:
                    pipeline.BATCH_NT = old

    def test_worker_error_propagates_and_no_partial_output(self):
        with self.assertRaises(CoalError):
            run(FIX / "codon_mismatch", "pep.fasta", ["nuc.fasta"], jobs=2, min_nt_per_job=0)
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out.fa"
            with self.assertRaises(CoalError):
                pipeline.run(FIX / "codon_mismatch" / "pep.fasta", [FIX / "codon_mismatch" / "nuc.fasta"],
                             out_path=str(out), jobs=2, min_nt_per_job=0)
            self.assertFalse(out.exists())
            self.assertEqual(list(Path(d).iterdir()), [])  # temp file cleaned up


if __name__ == "__main__":
    unittest.main()
