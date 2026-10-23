"""Reproductions of the findings of docs/stress-test-2026-10-23.md, as expected failures: each turns into an
"unexpected success" (remove the decorator) once the behaviour is fixed. Fast (< 2 s) unless COAL_SLOW is set."""
from __future__ import annotations

import os
import random
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CPP = os.environ.get("COAL_CPP_BIN") or str(ROOT / "build" / "coal-cpp")
HAVE_CPP = os.path.isfile(CPP) or os.path.isfile(CPP + ".exe")
ENV = dict(os.environ, COAL_ENGINE="python", PYTHONPATH=str(ROOT / "src"))


def run(engine, args, cwd):
    if engine == "python":
        cmd = [sys.executable, "-m", "coal", *args]
    else:
        cmd = [CPP, *args]
    return subprocess.run(cmd, cwd=cwd, capture_output=True, env=ENV, timeout=600)


def mixed_endings_case(d: Path, first_record_bytes: int):
    """Two genes; the DNA's first record is `first_record_bytes` long with LF line ends, the second record uses lone
    CR (classic Mac) line ends. PAL2NAL v14 reads this file; COAL decides the line-end style from the first 64 kB."""
    rng = random.Random(1)
    g1, g2 = "ATG" + "GCT" * 9, "ATG" + "TGG" * 9
    utr = "".join(rng.choice("ACGT") for _ in range(first_record_bytes))
    (d / "pep.fasta").write_bytes(b">a\nMAAAAAAAAA\n>b\nMWWWWWWWWW\n")
    (d / "nuc.fasta").write_bytes(b">a\n" + (g1 + utr).encode() + b"\n>b\r" + g2.encode() + b"\r")
    return b">a\n" + g1.encode() + b"\n>b\n" + g2.encode() + b"\n"


class MixedLineEndings(unittest.TestCase):
    """Stress-test finding 2: a lone-CR record after the first 64 kB of a DNA file is not seen as a record."""

    def check(self, engine, first_record_bytes):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            expected = mixed_endings_case(d, first_record_bytes)
            p = run(engine, ["pep.fasta", "nuc.fasta", "-output", "fasta"], tmp)
            self.assertEqual((p.returncode, p.stdout), (0, expected), p.stderr[:300])

    def test_lone_cr_record_within_64kb_python(self):
        self.check("python", 1000)

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_lone_cr_record_within_64kb_cpp(self):
        self.check("cpp", 1000)

    @unittest.expectedFailure
    def test_lone_cr_record_after_64kb_python(self):
        self.check("python", 70_000)

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    @unittest.expectedFailure
    def test_lone_cr_record_after_64kb_cpp(self):
        self.check("cpp", 70_000)


def repeat_case(d: Path, n: int):
    """A poly-A gene of n residues ending in W, placed after n more codons of poly-A in its DNA (an exact fit, found
    after a long low-complexity stretch)."""
    rng = random.Random(2)
    ala = ["GCA", "GCC", "GCG", "GCT"]
    (d / "pep.fasta").write_text(">r\n" + "A" * (n - 1) + "W\n")
    (d / "nuc.fasta").write_text(">r\n" + "".join(rng.choice(ala) for _ in range(n)) + "".join(rng.choice(ala) for _ in range(n - 1)) + "TGG\n")


@unittest.skipUnless(os.environ.get("COAL_SLOW"), "set COAL_SLOW=1 (about 10 s)")
@unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
class QuadraticExactFitScan(unittest.TestCase):
    """Stress-test finding 1: coal-cpp's exact-fit search (core.cpp find) tries every start and verifies each one
    codon by codon, so a low-complexity gene after a low-complexity stretch costs time quadratic in its length
    (4x the length: about 16x the time; the Python engine, which matches translated frames, grows linearly)."""

    @unittest.expectedFailure
    def test_time_grows_linearly_with_gene_length(self):
        times = {}
        for n in (10_000, 20_000):
            with tempfile.TemporaryDirectory() as tmp:
                repeat_case(Path(tmp), n)
                t = time.monotonic()
                p = run("cpp", ["pep.fasta", "nuc.fasta", "-o", "out.fa"], tmp)
                times[n] = time.monotonic() - t
                self.assertEqual(p.returncode, 0, p.stderr[:300])
        ratio = times[20_000] / max(times[10_000], 1e-3)
        self.assertLess(ratio, 3.0, f"2x the gene length took {ratio:.1f}x the time ({times})")


if __name__ == "__main__":
    unittest.main()
