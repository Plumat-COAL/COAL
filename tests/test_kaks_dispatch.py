"""YN / MYN: the AVX2 copy of the transition-probability matrix gives the same bits as the baseline copy.

On x86 with GCC or Clang, cpp/src/kaks.cpp compiles pmat_codon (the 61 x 61 eigendecomposition of each round of
YN and MYN, most of their time) twice, for the baseline ISA (SSE2) and for AVX2, and uses the AVX2 copy when the CPU
has it (a run-time check; no build or package needs an -m flag). Both copies do the same floating-point operations
in the same order, so the outputs must be byte-identical. `--kaks-baseline` (not in the help; testing only) makes a
run use the baseline copy: this test compares the two at 17 significant digits, both output formats, both YN
flavours and weightings, several genetic codes. Where the CPU has no AVX2 (or the build has no AVX2 copy: MSVC,
ARM) both runs take the baseline and the comparison holds trivially. Needs build/coal-cpp (COAL_CPP_BIN).
"""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import needs_cpp, run_cpp  # noqa: E402
from test_kaks_yn_shortcut import low_divergence_pairs, multi_difference  # noqa: E402

KAKS = ROOT / "tests" / "kaks"


@needs_cpp
class Dispatch(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def both(self, sub, path, *args):
        base = ["-i", path] if sub == "kaks" else [path]
        a = run_cpp(sub, *base, "-q", "--digits", "17", "-j", "2", *args)
        b = run_cpp(sub, *base, "-q", "--digits", "17", "-j", "2", "--kaks-baseline", *args)
        self.assertEqual(a.returncode, 0, a.stderr)
        self.assertEqual(b.returncode, 0, b.stderr)
        return a.stdout, b.stdout

    def test_same_bytes_on_both_paths(self):
        pairs = [p for p in low_divergence_pairs(600, 23) if multi_difference(*p)]   # (the matrix is computed)
        self.assertGreater(len(pairs), 120)
        gen = self.tmp / "multi.axt"
        gen.write_text("".join(f"p{k}&q{k}\n{a}\n{b}\n\n" for k, (a, b) in enumerate(pairs)), encoding="utf-8")
        for path in (gen, KAKS / "comeron_li_cases.axt", KAKS / "edge_cases.axt", KAKS / "sim_small.axt"):
            for flavour in ("kaks", "paml"):
                for weighting in ("1", "0"):
                    for sub in ("dnds", "kaks"):
                        with self.subTest(file=path.name, flavour=flavour, weighting=weighting, sub=sub):
                            fast, base = self.both(sub, path, "-m", "YN", "-m", "MYN", "--yn-flavour", flavour,
                                                   "--yn-weighting", weighting)
                            self.assertEqual(fast, base)
                            self.assertGreater(fast.count(b"\n"), 10)
        for code in ("2", "11"):
            with self.subTest(code=code):
                fast, base = self.both("dnds", gen, "-m", "YN", "-m", "MYN", "-c", code)
                self.assertEqual(fast, base)

    def test_every_method_on_both_paths(self):
        fast, base = self.both("dnds", KAKS / "sim_small.axt", *[x for m in ("YN", "NG", "LWL", "MLWL", "LPB", "MLPB",
                                                                             "MYN", "Comeron", "Li") for x in ("-m", m)])
        self.assertEqual(fast, base)

    def test_the_option_is_not_advertised(self):
        self.assertNotIn(b"kaks-baseline", run_cpp("dnds", "--help").stdout)
        self.assertNotIn(b"kaks-baseline", run_cpp("kaks", "--help").stdout)


if __name__ == "__main__":
    unittest.main()
