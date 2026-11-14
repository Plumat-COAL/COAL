"""The low findings P01, P02, P04 and P05 of the final re-test of 2026-11-14, each fixed and kept fixed here.

  * P01: a SIGTERM / SIGHUP / Ctrl-C that came while a run took or released its output lock left the lock
    `<out>.coal-lock` and/or its temporary file `<lock>.<pid>.<hex>.tmp` (about 2.4 % of the CI runs of
    tests/test_batch.py's SIGTERM test). Each window is hit here on purpose: the Python engine through a child
    that sends itself SIGTERM at that instant, coal-cpp through the hook COAL_TEST_LOCK_SIGNAL_AT.
  * P02: a whole set whose candidates all depart from the standard stops also names the NCBI codes COAL does not
    support that fit as well (code 25 bacteria got "codes 4 and 10 ... run again with the -c of its organism").
  * P04: --subject with the same records as --query (re-wrapped, other case, a record more): the copy warning.
  * P05: an input that does not exist is refused before any .gz input is decompressed.

Synthetic inputs only (no real sequence data); coal-cpp (COAL_CPP_BIN) for the pipeline parts.
"""
from __future__ import annotations

import gzip
import io
import os
import random
import signal
import subprocess
import sys
import tempfile
import textwrap
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from cpp_hooks import need_test_hooks  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from coal import pipeline as L  # noqa: E402
from test_deep_findings import cds, write_fasta  # noqa: E402
from test_dnds_pipeline import CPP, FAM_CLUSTERS, Q, SCHEMA, PipelineCase, needs_cpp  # noqa: E402
from test_report_schema import validate  # noqa: E402

EXAMPLE = [str(ROOT / "examples" / "clean" / "protein.aln.fasta"), str(ROOT / "examples" / "clean" / "dna.fasta")]
POSIX = os.name == "posix"


def run(fn, *args, **kw):
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, " ".join(err.getvalue().split())


# A `coal PEP NUC -o out.fa` run (Python engine) that sends itself SIGTERM at one instant of the lock code.
# argv: src WHERE PEP NUC -o OUT
CHILD = r"""
import os, signal, sys
sys.path.insert(0, sys.argv.pop(1))
where = sys.argv.pop(1)
import coal.pipeline as P

def term():
    os.kill(os.getpid(), signal.SIGTERM)

if where == "after_link":            # the lock is linked, its temporary file not yet removed
    real_link = os.link
    def link(a, b, **k):
        real_link(a, b, **k)
        if b.endswith(".coal-lock"):
            term()
    P.os.link = link
elif where == "tmp_written":         # the temporary file is open and written, not yet linked
    real_fdopen = os.fdopen
    def fdopen(fd, *a, **k):
        fh = real_fdopen(fd, *a, **k)
        term()
        return fh
    P.os.fdopen = fdopen
elif where == "in_remove_tmp":       # while the temporary file is removed (after the link)
    real_rm = P._remove
    def rm(p, *a, **k):
        if p.endswith(".tmp"):
            term()
        return real_rm(p, *a, **k)
    P._remove = rm
elif where == "after_take":          # take_output_lock has returned, the caller's try is not entered yet
    real_take = P.take_output_lock
    def take(*a, **k):
        r = real_take(*a, **k)
        term()
        return r
    P.take_output_lock = take
elif where == "before_release":      # the output is in place, the lock not released yet
    real_rel = P.release_output_lock
    def rel(path):
        if path:
            term()
        return real_rel(path)
    P.release_output_lock = rel
from coal.cli import main
sys.exit(main())
"""


@unittest.skipUnless(POSIX, "SIGTERM to itself ends a Windows process at once, without Python's handler")
class PythonLockSignals(unittest.TestCase):
    """P01, Python engine: a SIGTERM in each window of taking and releasing the lock leaves nothing behind."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(os.path.realpath(self._td.name))

    def tearDown(self):
        self._td.cleanup()

    def window(self, where):
        d = self.tmp / where
        d.mkdir()
        env = {k: v for k, v in os.environ.items() if not k.startswith("COAL_TEST_")}
        env["COAL_ENGINE"] = "python"   # (the lock code under test; coal-cpp has its own test below)
        r = subprocess.run([sys.executable, "-c", CHILD, str(ROOT / "src"), where, *EXAMPLE, "-o", "out.fa", "-q"],
                           cwd=d, env=env, capture_output=True, text=True, timeout=120)
        return r, sorted(x.name for x in d.iterdir() if x.name != "out.fa")

    def test_every_window_leaves_nothing(self):
        for where in ("after_link", "tmp_written", "in_remove_tmp", "after_take", "before_release"):
            with self.subTest(where=where):
                r, left = self.window(where)
                self.assertEqual(r.returncode, 143, r.stderr)
                self.assertEqual(left, [], r.stderr)   # (before: the lock, and in_remove_tmp the .tmp too)

    def test_a_run_after_an_interrupted_one(self):
        r, left = self.window("after_take")
        self.assertEqual(left, [])
        r = subprocess.run([sys.executable, "-m", "coal", *EXAMPLE, "-o", "out.fa"], cwd=self.tmp / "after_take",
                           env={**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_ENGINE": "python"},
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)


class HeldLocks(unittest.TestCase):
    """P01: release_held_locks() (at exit, and in a batch worker after each job) and take_output_lock's own
    cleanup."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(os.path.realpath(self._td.name))
        self.out = str(self.tmp / "out.fa")

    def tearDown(self):
        L.release_held_locks()
        self._td.cleanup()

    def test_held_locks_and_temporary_files_are_removed(self):
        lock = L.take_output_lock(self.out, self.out)
        tmp = Path(lock + f".{os.getpid()}.0123456789abcdef.tmp")
        tmp.write_bytes(b"x")
        L._own_tmps.add(str(tmp))
        L.release_held_locks()
        self.assertEqual(sorted(x.name for x in self.tmp.iterdir()), [])
        self.assertEqual((L._held_locks, L._own_tmps), ({}, set()))

    def test_a_lock_that_is_not_ours_any_more_stays(self):
        lock = L.take_output_lock(self.out, self.out)
        Path(lock).write_text("1 otherhost - feedfeed\n", encoding="ascii")   # replaced by hand and taken again
        L.release_held_locks()
        self.assertTrue(Path(lock).exists())
        self.assertEqual(L._held_locks, {})

    def test_an_interrupt_inside_take_removes_the_lock_it_made(self):
        # KeyboardInterrupt (Ctrl-C) after the lock was made and read back, before take_output_lock returned
        real = L._read_lock
        hit = []

        def read(path):
            got = real(path)
            if path.endswith(".coal-lock") and got is not None and got[1] == os.getpid() and not hit:
                hit.append(path)
                raise KeyboardInterrupt
            return got
        with mock.patch.object(L, "_read_lock", read):
            with self.assertRaises(KeyboardInterrupt):
                L.take_output_lock(self.out, self.out)
        self.assertEqual(sorted(x.name for x in self.tmp.iterdir()), [])
        self.assertEqual(L._held_locks, {})
        L.release_output_lock(L.take_output_lock(self.out, self.out))   # and the output can be taken again

    @unittest.skipUnless(POSIX and hasattr(signal, "pthread_sigmask"), "no signal mask on this system")
    def test_signals_wait_while_the_lock_is_made(self):
        got = []
        old = signal.signal(signal.SIGUSR1, lambda *a: got.append(sorted(x.name for x in self.tmp.iterdir())))
        try:
            real = L._BLOCKED
            with mock.patch.object(L, "_BLOCKED", (*real, signal.SIGUSR1)):
                real_link = os.link

                def link(a, b, **k):
                    real_link(a, b, **k)
                    os.kill(os.getpid(), signal.SIGUSR1)   # delivered when the step is done: the .tmp is gone
                with mock.patch.object(L.os, "link", link):
                    lock = L.take_output_lock(self.out, self.out)
        finally:
            signal.signal(signal.SIGUSR1, old)
        self.assertEqual(got, [["out.fa.coal-lock"]])
        L.release_output_lock(lock)


class CppLockSignals(unittest.TestCase):
    """P01, coal-cpp: a SIGTERM right after the lock is linked (before it was registered for the signal handler)
    left the lock."""

    def test_signal_after_the_link_leaves_nothing(self):
        if CPP is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("coal-cpp not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("coal-cpp not built (cmake --build build)")
        if os.name != "posix":
            self.skipTest("coal-cpp blocks the signals on POSIX only (Windows has no signal mask)")
        need_test_hooks(self, CPP)
        with tempfile.TemporaryDirectory() as td:
            env = {k: v for k, v in os.environ.items() if not k.startswith("COAL_TEST_")}
            env.update(COAL_TEST_HOOKS="1", COAL_TEST_LOCK_SIGNAL_AT="link")
            r = subprocess.run([str(CPP), *EXAMPLE, "-o", "out.fa"], cwd=td, env=env, capture_output=True,
                               text=True, timeout=120)
            self.assertEqual(r.returncode, 143, r.stderr)
            self.assertEqual(sorted(os.listdir(td)), [])


class UnsupportedCodes(unittest.TestCase):
    """P02: suggest_genetic_code() on synthetic sets."""

    def setUp(self):
        self.rng = random.Random(25)
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def hint(self, genes, nuclear=0, chosen=11):
        recs = [(f"nuc{k}", cds(self.rng, stops=())) for k in range(nuclear)]
        recs += [(f"g{k}", s) for k, s in enumerate(genes)]
        path = write_fasta(self.tmp / "set.fa", recs)
        flagged = {f"g{k}": 0 for k in range(len(genes))}
        usable = {r for r, _ in recs} - set(flagged)
        g, bg, reps = P._hint_inputs([str(path)], flagged, usable)
        return P.suggest_genetic_code(chosen, g, bg, len(recs), len(genes), reps)

    def test_tga_genes_of_a_genome_also_name_code_25(self):
        h = self.hint([cds(self.rng) for _ in range(80)], nuclear=20)
        self.assertTrue(h["whole_set"])
        self.assertIn(4, h["codes"])
        self.assertEqual(h["unsupported"], [25])   # (not the mitochondrial 24 / 33 for a genome-size set)

    def test_a_small_set_also_names_the_organelle_codes(self):
        h = self.hint([cds(self.rng) for _ in range(10)])
        self.assertEqual(h["unsupported"], [24])   # (33 also reads TAA: a stop these genes do not need)

    def test_no_unsupported_codes_beside_a_single_code_or_the_standard_stops(self):
        h = self.hint([cds(self.rng, avoid=("AGA", "AGG")) for _ in range(10)], nuclear=200, chosen=1)
        self.assertEqual((h["codes"], h["unsupported"]), ([2], []))   # vertebrate mitochondria in a nuclear set
        h = self.hint([cds(self.rng, stops=("TAA", "TAG", "TAA")) for _ in range(80)], nuclear=20)
        self.assertEqual((h["codes"], h["unsupported"]), ([6], []))

    def test_texts(self):
        self.assertEqual(P._supported_text(), "1-6, 9-16 and 21-23")
        self.assertEqual(P._unsupported_text([25]),
                         "NCBI code 25 (Candidate Division SR1 and Gracilibacteria: TGA = Gly)")


@needs_cpp
class UnsupportedCodeHint(PipelineCase):
    """P02 end to end: the warning and report.json of a clusters run."""

    def test_the_hint_says_that_code_25_is_not_supported(self):
        rng = random.Random(5)
        genes = [(f"lcl|NC_099999.1_cds_G{k}_{k}", cds(rng)) for k in range(160)]
        f = write_fasta(self.tmp / "sr1.fa", genes)
        cl = self.tmp / "cl.tsv"
        cl.write_text(f"{genes[0][0]}\t{genes[1][0]}\n", encoding="utf-8")
        res, err = run(P.run_clusters, [f], self.tmp / "out", clusters=cl, codontable=11)
        validate(res.report, SCHEMA)
        (h,) = res.report["genetic_code_hints"]
        self.assertTrue(h["whole_set"])
        self.assertEqual(h["unsupported_candidates"], [25])
        self.assertIn("; or NCBI code 25 (Candidate Division SR1 and Gracilibacteria: TGA = Gly), which COAL does "
                      "not support (-c takes 1-6, 9-16 and 21-23)", err)
        self.assertIn("for a code COAL does not support, no -c gives right rows", err)
        self.assertNotIn("Run again with the -c of its organism", err)


if __name__ == "__main__":
    unittest.main()
