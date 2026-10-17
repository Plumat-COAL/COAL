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


if __name__ == "__main__":
    unittest.main()
