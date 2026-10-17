"""Ordered streaming thread pool of the C++ engine (CPP-5).

Through `coal-cpp-dump pool N JOBS BATCH [FAIL_AT]` (cpp/src/pool.hpp with a stand-in per-item
function whose items finish out of order) and through coal-cpp itself:
  * output always in input order, the same with -j 1, 2, 4 and the default (hardware threads);
  * at most two batches in memory (the batch being worked on + the one being read);
  * -j N uses N threads per batch (fewer if the batch is smaller); default = hardware threads;
  * a failure in the per-item function is reported on the calling thread (exit 1, no crash) and
    nothing after the last complete batch is written;
  * coal-cpp gives identical output for every -j with small batches (COAL_BATCH_NT) on 300 sequences.
Skipped if the binaries are not built.
"""
from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_cpp_readers import ROOT, TIMEOUT, cpp_dump  # noqa: E402
from test_cpp_cli import cpp_binary  # noqa: E402

M64 = (1 << 64) - 1


def stand_in(i: int) -> int:
    """The dump's stand-in function for an item that is not first in its batch."""
    x = i
    for _ in range((i % 50 + 1) * 200):
        x = (x * 6364136223846793005 + 1442695040888963407) & M64
    return x


def pool(dump, n, jobs, batch, fail_at=None):
    args = [str(dump), "pool", str(n), str(jobs), str(batch)] + ([str(fail_at)] if fail_at is not None else [])
    r = subprocess.run(args, capture_output=True, text=True, timeout=TIMEOUT)
    lines = r.stdout.splitlines()
    stats = json.loads(lines[-1]) if lines else {}
    items = [tuple(map(int, ln.split())) for ln in lines[:-1]]
    return r.returncode, items, stats


class PoolTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dump = cpp_dump()
        if cls.dump is None:
            raise unittest.SkipTest("coal-cpp-dump not built")

    def test_order_and_results_independent_of_jobs(self):
        n, batch = 2000, 100
        rc, ref, st = pool(self.dump, n, 1, batch)
        self.assertEqual(rc, 0)
        self.assertEqual([i for i, _ in ref], list(range(n)))
        self.assertFalse(st["reordered"])
        for i in (1, 7, 199, 1234, 1999):  # spot check against Python
            self.assertEqual(ref[i][1], stand_in(i))
        for jobs in (2, 4, 0):
            with self.subTest(jobs=jobs):
                rc, got, st = pool(self.dump, n, jobs, batch)
                self.assertEqual(rc, 0)
                self.assertEqual(got, ref)
                self.assertEqual((st["items"], st["batches"]), (n, n // batch))

    def test_items_finish_out_of_order_but_are_written_in_order(self):
        rc, got, st = pool(self.dump, 400, 4, 50)
        self.assertEqual(rc, 0)
        self.assertTrue(st["reordered"])  # completion order differed from input order...
        self.assertEqual([i for i, _ in got], list(range(400)))  # ...output did not

    def test_memory_bounded_by_two_batches(self):
        for n, batch in ((1000, 10), (1000, 64), (5, 64), (130, 64)):
            with self.subTest(n=n, batch=batch):
                rc, got, st = pool(self.dump, n, 4, batch)
                self.assertEqual(rc, 0)
                self.assertEqual(len(got), n)
                self.assertLessEqual(st["max_live"], 2 * batch)
                self.assertEqual(st["batches"], -(-n // batch))

    def test_threads(self):
        for jobs, batch, want in ((1, 50, 1), (2, 50, 2), (4, 50, 4), (8, 3, 3)):
            with self.subTest(jobs=jobs, batch=batch):
                self.assertEqual(pool(self.dump, 200, jobs, batch)[2]["max_threads"], want)
        st = pool(self.dump, 200, 0, 50)[2]
        self.assertEqual(st["max_threads"], min(st["threads_default"], 50))
        if os.cpu_count():
            self.assertEqual(st["threads_default"], os.cpu_count())

    def test_empty_and_single(self):
        self.assertEqual(pool(self.dump, 0, 4, 10)[1:], ([], {"batches": 0, "items": 0, "max_live": 0,
                                                              "max_threads": 0, "threads_default":
                                                                  pool(self.dump, 0, 4, 10)[2]["threads_default"],
                                                              "reordered": False}))
        rc, got, st = pool(self.dump, 1, 4, 10)
        self.assertEqual((rc, [i for i, _ in got], st["max_threads"]), (0, [0], 1))

    def test_failure_is_reported_not_crashed(self):
        for jobs in (1, 4):
            with self.subTest(jobs=jobs):
                rc, items, st = pool(self.dump, 10, jobs, 4, fail_at=6)
                self.assertEqual(rc, 1)
                self.assertEqual(st, {"error": "stand-in failure at item 6"})
                self.assertEqual([i for i, _ in items], [0, 1, 2, 3])  # only the batch before the failure


AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
         "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


class EngineJobsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bin = cpp_binary()
        if cls.bin is None:
            raise unittest.SkipTest("coal-cpp not built")

    def test_same_output_for_every_jobs_value(self):
        rnd = random.Random(5)
        width = 120
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            pep, nuc = [], []
            for k in range(300):
                prot = "M" + "".join(rnd.choice(AA) for _ in range(rnd.randint(20, 90)))
                row = list(prot)
                while len(row) < width:
                    row.insert(rnd.randint(1, len(row)), "-")
                pep.append(f">s{k}\n{''.join(row)}\n")
                dna = "".join(CODON[a] for a in prot) + "TAA"
                if k % 7 == 0:
                    dna = dna[:30] + "TGG" + dna[33:]  # one wrong codon: a mismatch warning
                nuc.append(f">s{k}\n{dna}\n")
            (t / "pep.fa").write_text("".join(pep))
            (t / "nuc.fa").write_text("".join(nuc))
            env = dict(os.environ, COAL_BATCH_NT="500")  # ~3 sequences per batch: ~100 batches
            for fmt in ("fasta", "clustal", "paml", "codon"):
                outs = {}
                for jobs in ("1", "2", "4", "0"):
                    r = subprocess.run([str(self.bin), "pep.fa", "nuc.fa", "-output", fmt, "-j", jobs], cwd=t,
                                       env=env, capture_output=True, text=True, timeout=TIMEOUT)
                    outs[jobs] = (r.returncode, r.stdout, r.stderr)
                with self.subTest(fmt=fmt):
                    self.assertTrue(outs["1"][1])
                    self.assertIn("warning", outs["1"][2])  # the wrong codons were seen
                    for jobs in ("2", "4", "0"):
                        self.assertEqual(outs[jobs], outs["1"], jobs)
                    if fmt == "fasta":  # input order, across ~100 batches
                        ids = [ln[1:] for ln in outs["1"][1].splitlines() if ln.startswith(">")]
                        self.assertEqual(ids, [f"s{k}" for k in range(300)])


if __name__ == "__main__":
    unittest.main()
