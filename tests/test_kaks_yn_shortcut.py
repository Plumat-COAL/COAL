"""YN / MYN without the P(t) matrix where nothing reads it: the same bits as with it.

The transition-probability matrix of the codon model (a 61 x 61 eigendecomposition per round) only weights the
pathways between codons that differ at two or three positions. For a pair of sequences without such a codon it is not
computed (cpp/src/kaks.cpp needs_pathway_weights). `--yn-always-pmat` (not in the help; testing only) computes it all
the same, as every release before did: the two outputs must be byte-identical at 17 significant digits, in both
output formats, for both flavours, on every platform this runs on. docs/kaks.md "Benchmark" has the measurement on
the validation sets. Needs build/coal-cpp (COAL_CPP_BIN).
"""
from __future__ import annotations

import random
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import needs_cpp, run_cpp  # noqa: E402

KAKS = ROOT / "tests" / "kaks"
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SENSE = [c for c in CODONS if c not in ("TAA", "TAG", "TGA")]


def low_divergence_pairs(n: int, seed: int) -> list[tuple[str, str]]:
    """Pairs of 20 to 200 codons at several divergences, from identical to far apart: many of the close ones have
    no codon that differs at two positions."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        a = [rng.choice(SENSE) for _ in range(rng.randint(20, 200))]
        rate = rng.choice([0.0, 0.003, 0.01, 0.03, 0.1, 0.3])
        b = []
        for c in a:
            d = c
            for pos in range(3):
                if rng.random() < rate:
                    e = d[:pos] + rng.choice(BASES) + d[pos + 1:]
                    if e in SENSE:
                        d = e
            b.append(d)
        if rng.random() < 0.1:   # a gap codon and an ambiguous one
            b[len(b) // 2] = "---"
            a[len(a) // 3] = "NNN"
        out.append(("".join(a), "".join(b)))
    return out


def multi_difference(a: str, b: str) -> bool:
    for k in range(0, len(a) - 2, 3):
        x, y = a[k:k + 3], b[k:k + 3]
        if "-" in x + y or "N" in x + y:
            continue
        if sum(p != q for p, q in zip(x, y)) >= 2:
            return True
    return False


@needs_cpp
class Shortcut(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def both(self, sub, path, *args):
        base = ["-i", path] if sub == "kaks" else [path]
        a = run_cpp(sub, *base, "-q", "--digits", "17", "-j", "2", *args)
        b = run_cpp(sub, *base, "-q", "--digits", "17", "-j", "2", "--yn-always-pmat", *args)
        self.assertEqual(a.returncode, 0, a.stderr)
        self.assertEqual(b.returncode, 0, b.stderr)
        return a.stdout, b.stdout

    def test_same_bytes_with_and_without_the_matrix(self):
        pairs = low_divergence_pairs(300, 11)
        with_multi = sum(multi_difference(a, b) for a, b in pairs)
        self.assertGreater(with_multi, 60)                     # both kinds of pair are in the set
        self.assertGreater(len(pairs) - with_multi, 60)
        gen = self.tmp / "low.axt"
        gen.write_text("".join(f"p{k}&q{k}\n{a}\n{b}\n\n" for k, (a, b) in enumerate(pairs)), encoding="utf-8")
        for path in (gen, KAKS / "comeron_li_cases.axt", KAKS / "edge_cases.axt", KAKS / "sim_small.axt"):
            for flavour in ("kaks", "paml"):
                for sub in ("dnds", "kaks"):
                    with self.subTest(file=path.name, flavour=flavour, sub=sub):
                        fast, full = self.both(sub, path, "-m", "YN", "-m", "MYN", "--yn-flavour", flavour)
                        self.assertEqual(fast, full)
                        self.assertGreater(fast.count(b"\n"), 10)
        # ... and the other genetic codes take the same route
        fast, full = self.both("dnds", gen, "-m", "YN", "-m", "MYN", "-c", "2")
        self.assertEqual(fast, full)

    def test_iteration_count_and_statuses_are_those_of_the_full_computation(self):
        """A pair without a multi-difference codon still runs its two rounds: nothing about convergence changes."""
        p = self.tmp / "one.axt"
        p.write_text("a&b\nATGAAAGTTCTGGCAGGTCCGTTAGCTGCA\nATGAAGGTTCTGGCAGGACCGTTAGCTGCG\n\n", encoding="utf-8")
        fast, full = self.both("dnds", p, "-m", "YN", "-m", "MYN")
        self.assertEqual(fast, full)
        rows = [ln.split("\t") for ln in fast.decode().splitlines()[1:]]
        self.assertEqual([r[13] for r in rows], ["ok", "ok"])

    def test_the_option_is_not_advertised(self):
        self.assertNotIn(b"always-pmat", run_cpp("dnds", "--help").stdout)
        self.assertNotIn(b"always-pmat", run_cpp("kaks", "--help").stdout)


if __name__ == "__main__":
    unittest.main()
