"""Regression tests for the independent engine review (overnight/engine-fixes/CHANGES.md, ENG-a .. ENG-i).

Every input is built inside a temporary directory created by the test; -o is only
ever pointed at paths inside that directory.
"""
from __future__ import annotations

import os
import random
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
TIMEOUT = 60

PEP = ">s1\nMKV-AL\n>s2\nMKVSAL\n"
NUC = ">s1\nATGAAAGTTGCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n"
EXPECTED = ">s1\nATGAAAGTT---GCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n"


def coal(args, cwd, stdin=None, timeout=TIMEOUT, **kw):
    env = dict(os.environ, PYTHONPATH=str(SRC))
    return subprocess.run([sys.executable, "-m", "coal", *args], cwd=cwd, env=env, input=stdin,
                          capture_output=True, text=True, timeout=timeout, **kw)


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_review_")
        self.d = Path(self._td.name)
        (self.d / "p.fa").write_text(PEP)
        (self.d / "n.fa").write_text(NUC)

    def tearDown(self):
        self._td.cleanup()

    def leftovers(self):
        return sorted(p.name for p in self.d.iterdir() if p.name.startswith(".coal_"))


# ---------------------------------------------------------------- ENG-a
class TestEngA_OutputPath(TempDirCase):
    def test_symlink_is_written_through(self):
        target = self.d / "target.fa"
        target.write_text("old\n")
        link = self.d / "link.fa"
        link.symlink_to(target)
        r = coal(["p.fa", "n.fa", "-o", "link.fa"], self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(link.is_symlink())
        self.assertEqual(target.read_text(), EXPECTED)
        self.assertEqual(self.leftovers(), [])

    def test_new_file_mode_follows_umask(self):
        old = os.umask(0o022)
        try:
            r = coal(["p.fa", "n.fa", "-o", "out.fa"], self.d)
        finally:
            os.umask(old)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(stat.S_IMODE((self.d / "out.fa").stat().st_mode), 0o644)

    @unittest.skipUnless(hasattr(os, "mkfifo"), "needs FIFOs")
    def test_fifo_is_written_not_replaced(self):
        fifo = self.d / "out.fifo"
        os.mkfifo(fifo)
        reader = subprocess.Popen(["cat", str(fifo)], stdout=subprocess.PIPE, text=True)
        try:
            r = coal(["p.fa", "n.fa", "-o", str(fifo)], self.d)
            out, _ = reader.communicate(timeout=20)
        finally:
            if reader.poll() is None:
                reader.kill()
                reader.wait()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out, EXPECTED)
        self.assertTrue(stat.S_ISFIFO(fifo.stat().st_mode))
        self.assertEqual(self.leftovers(), [])


# ---------------------------------------------------------------- ENG-b
class TestEngB_StreamInputs(TempDirCase):
    def test_protein_from_pipe(self):
        r = coal(["/dev/stdin", "n.fa"], self.d, stdin=PEP)
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)

    def test_dna_from_pipe(self):
        tmp = Path(tempfile.gettempdir())
        before = set(tmp.glob(".coal_in_*"))
        r = coal(["p.fa", "/dev/stdin"], self.d, stdin=NUC)
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)
        self.assertEqual(set(tmp.glob(".coal_in_*")) - before, set())  # spooled copy removed

    @unittest.skipUnless(hasattr(os, "mkfifo"), "needs FIFOs")
    def test_fifo_inputs(self):  # what <(cat p.fa) gives
        os.mkfifo(self.d / "pf")
        os.mkfifo(self.d / "nf")
        writer = subprocess.Popen(["sh", "-c", "cat p.fa > pf & cat n.fa > nf & wait"], cwd=self.d)
        try:
            r = coal(["pf", "nf"], self.d)
        finally:
            try:
                writer.wait(timeout=5)
            except subprocess.TimeoutExpired:
                writer.kill()
                writer.wait()
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)


# ---------------------------------------------------------------- ENG-c
def _random_pair(seed=3, n=60):
    rng = random.Random(seed)
    prot = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n - 1))
    dna = "".join(rng.choice("ACGT") for _ in range(3 * n))
    return f">g1\n{prot}\n>g2\n{prot}\n", f">g1\n{dna}\n>g2\n{dna}\n"


class TestEngC_Unrelated(TempDirCase):
    def test_unrelated_protein_and_dna_is_an_error(self):
        pep, nuc = _random_pair()
        (self.d / "rp.fa").write_text(pep)
        (self.d / "rn.fa").write_text(nuc)
        for extra in ([], ["-q"]):
            r = coal(["rp.fa", "rn.fa", *extra], self.d)
            self.assertNotEqual(r.returncode, 0, r.stdout)
            self.assertEqual(r.stdout, "")
            self.assertIn("g1: protein and DNA are inconsistent", r.stderr)

    def test_no_anchor_matched_warns_even_with_q(self):
        # Golden table3/fasta: PAL2NAL succeeds although no 10-residue anchor matches
        # (7 of 25 codons differ under table 1), so this stays a success, with a notice.
        fx = ROOT / "tests" / "fixtures" / "table3"
        gold = (ROOT / "tests" / "golden" / "table3" / "fasta.stdout").read_text()
        r = coal([str(fx / "pep.fasta"), str(fx / "nuc.fasta"), "-q"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, gold), r.stderr)
        self.assertIn("s0: no 10-residue segment", r.stderr)
        self.assertIn("7 of 25", r.stderr)
        self.assertNotIn("pepAlnPos", r.stderr)


# ---------------------------------------------------------------- ENG-d
class TestEngD_PartialIdMatch(TempDirCase):
    def test_some_ids_match_warns_and_pairs_by_order(self):
        (self.d / "ip.fa").write_text(">geneA\nMKVAL\n>geneB\nMKVAL\n")
        (self.d / "in.fa").write_text(">geneB\nATGAAAGTTGCTCTG\n>geneX\nATGAAGGTCGCACTA\n")
        for extra in ([], ["-q"]):
            r = coal(["ip.fa", "in.fa", *extra], self.d)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(r.stdout, ">geneA\nATGAAAGTTGCTCTG\n>geneB\nATGAAGGTCGCACTA\n")
            self.assertIn("1 of 2 protein IDs", r.stderr)
            self.assertIn("'geneA' gets DNA 'geneB'", r.stderr)

    def test_no_ids_match_stays_silent(self):  # golden id_mismatch_same_order: PAL2NAL pairs by order
        (self.d / "ip.fa").write_text(">a\nMKVAL\n>b\nMKVAL\n")
        (self.d / "in.fa").write_text(">x\nATGAAAGTTGCTCTG\n>y\nATGAAGGTCGCACTA\n")
        r = coal(["ip.fa", "in.fa"], self.d)
        self.assertEqual((r.returncode, r.stderr), (0, ""))


# ---------------------------------------------------------------- ENG-e
class TestEngE_LongFallback(unittest.TestCase):
    def test_100k_codons_with_one_n_is_fast(self):
        sys.path.insert(0, str(SRC))
        import time
        from coal.core import thread_one
        from coal.tables import get_code
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*")
        rng = random.Random(5)
        cod = ["ATG"] + [rng.choice(sense) for _ in range(99_999)]
        prot = "".join(code[c] for c in cod)
        dna = "".join(cod)
        k = 3 * 50_000 + 1
        dna = dna[:k] + "N" + dna[k + 1:]
        t = time.perf_counter()
        r = thread_one("g", prot, dna, 1)
        dt = time.perf_counter() - t
        self.assertEqual(r.dna, dna)
        self.assertEqual(r.mismatch_cols, {50_001})
        self.assertLess(dt, 1.0)


# ---------------------------------------------------------------- ENG-f
class TestEngF_Interrupts(TempDirCase):
    """coal is stopped while it copies its finished temp file to a stdout pipe that the
    test does not drain (so it is blocked mid-copy, with the temp file present)."""

    def setUp(self):
        super().setUp()
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(200, 1000, self.d / "big", seed=3)   # ~600 kB of output, > a pipe buffer
        self.tmpdir = self.d / "tmp"
        self.tmpdir.mkdir()

    def start(self):
        env = dict(os.environ, PYTHONPATH=str(SRC), TMPDIR=str(self.tmpdir))
        p = subprocess.Popen([sys.executable, "-m", "coal", "big/pep.fasta", "big/nuc.fasta", "-j", "1"],
                             cwd=self.d, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertTrue(p.stdout.readline().startswith(b">seq0"))  # now in the copy phase
        return p

    def stop(self, p, sig):
        p.send_signal(sig)
        try:
            p.stdout.close()  # only after the signal, so no broken pipe first
            err = p.stderr.read().decode()
            p.wait(timeout=TIMEOUT)
        finally:
            if p.poll() is None:
                p.kill()
                p.wait()
        return p.returncode, err

    def test_sigterm_cleans_up(self):
        import signal
        p = self.start()
        rc, err = self.stop(p, signal.SIGTERM)
        self.assertEqual(rc, 128 + signal.SIGTERM, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(list(self.tmpdir.iterdir()), [])

    def test_sigint_message_and_130(self):
        import signal
        p = self.start()
        rc, err = self.stop(p, signal.SIGINT)
        self.assertEqual(rc, 130, err)
        self.assertNotIn("Traceback", err)
        self.assertIn("coal: interrupted", err)
        self.assertEqual(list(self.tmpdir.iterdir()), [])

    def test_broken_pipe_is_quiet(self):  # coal ... | head -1
        from gen_data import generate
        generate(2000, 1000, self.d / "huge", seed=3)   # ~6 MB of output, > any pipe buffer
        env = dict(os.environ, PYTHONPATH=str(SRC), TMPDIR=str(self.tmpdir))
        p = subprocess.Popen([sys.executable, "-m", "coal", "huge/pep.fasta", "huge/nuc.fasta", "-j", "1"],
                             cwd=self.d, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            self.assertTrue(p.stdout.readline().startswith(b">seq0"))
            p.stdout.close()
            err = p.stderr.read().decode()
            p.wait(timeout=TIMEOUT)
        finally:
            if p.poll() is None:
                p.kill()
                p.wait()
        self.assertEqual(err, "")
        self.assertEqual(p.returncode, 141)
        self.assertEqual(list(self.tmpdir.iterdir()), [])


# ---------------------------------------------------------------- ENG-g
class TestEngG_Table3Starts(TempDirCase):
    def test_gtg_start_in_table3(self):
        (self.d / "gp.fa").write_text(">g\nMKVAS\n")
        (self.d / "gn.fa").write_text(">g\nGTGAAAGTTGCTTCT\n")
        r = coal(["gp.fa", "gn.fa", "-codontable", "3"], self.d)
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, ">g\nGTGAAAGTTGCTTCT\n", ""))

    def test_table3_start_set(self):
        sys.path.insert(0, str(SRC))
        from coal.tables import get_code
        self.assertEqual(set(get_code(3).starts), {"ATA", "ATG", "GTG"})


if __name__ == "__main__":
    unittest.main()
