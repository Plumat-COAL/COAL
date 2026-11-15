"""YN / MYN: the AVX2 copy of the transition-probability matrix gives the same bits as the baseline copy.

On x86 with GCC or Clang, cpp/src/kaks.cpp compiles pmat_codon (the 61 x 61 eigendecomposition of each round of
YN and MYN, most of their time) twice, for the baseline ISA (SSE2) and for AVX2, and uses the AVX2 copy when the CPU
has it (a run-time check; no build or package needs an -m flag). Both copies do the same floating-point operations
in the same order, so the outputs must be byte-identical. `--kaks-baseline` (not in the help; testing only) makes a
run use the baseline copy, and `--kaks-isa` (testing only too) prints which copy a run would use. This test first
checks that the AVX2 copy is the one used (where the CPU has AVX2 and the build has the copy), then compares the two
at 17 significant digits: both output formats, both YN flavours and weightings, several genetic codes. Where there
is no AVX2 copy to compare (a CPU without AVX2; ARM, where it is not compiled; MSVC, which has no target attribute)
the comparison is skipped with that reason, so a pass always means both paths ran. Needs build/coal-cpp
(COAL_CPP_BIN).
"""
from __future__ import annotations

import os
import platform
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import needs_cpp, run_cpp  # noqa: E402
from test_kaks_yn_shortcut import low_divergence_pairs, multi_difference  # noqa: E402

KAKS = ROOT / "tests" / "kaks"


def cpu_has_avx2() -> bool | None:
    """True / False where the operating system says, None where this cannot tell."""
    if platform.machine().lower() not in ("x86_64", "amd64", "i386", "i686", "x86"):
        return False
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/cpuinfo", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith("flags"):
                        return "avx2" in line.split()
            return None
        if sys.platform == "darwin":
            r = subprocess.run(["sysctl", "-n", "machdep.cpu.leaf7_features"], capture_output=True, text=True)
            return "AVX2" in r.stdout.split() if r.returncode == 0 else None
    except OSError:
        return None
    return None   # (Windows: the CI build there is MSVC's, which has no AVX2 copy)


@needs_cpp
class Dispatch(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def isa(self, *args) -> str:
        r = run_cpp("dnds", "--kaks-isa", *args)
        self.assertEqual(r.returncode, 0, r.stderr)
        return r.stdout.decode().strip()

    def require_avx2_path(self):
        """Skip (with the reason) unless runs without --kaks-baseline take the AVX2 copy."""
        isa, avx2 = self.isa(), cpu_has_avx2()
        if isa == "avx2":
            self.assertNotEqual(avx2, False, "the AVX2 copy runs on a CPU said to have no AVX2")
            return
        self.assertEqual(isa, "baseline")
        if avx2 is False:
            self.skipTest("this CPU has no AVX2: only the baseline copy can run")
        if avx2 is None:
            self.skipTest("cannot tell whether this CPU has AVX2, and the build runs the baseline copy")
        if os.name == "nt":
            self.skipTest("this build has no AVX2 copy (MSVC has no target attribute; cpp/src/kaks.cpp COAL_KAKS_AVX2)")
        self.fail("the CPU has AVX2 but coal-cpp runs the baseline copy: the AVX2 copy is missing from this "
                  "GCC / Clang build on x86 (cpp/src/kaks.cpp COAL_KAKS_AVX2)")

    def test_which_copy_runs(self):
        self.assertEqual(self.isa("--kaks-baseline"), "baseline")
        self.assertEqual(self.isa("--kaks-isa", "--kaks-baseline"), "baseline")
        if platform.machine().lower() in ("arm64", "aarch64"):
            self.assertEqual(self.isa(), "baseline")   # (no AVX2 copy is compiled there)
        self.require_avx2_path()

    def both(self, sub, path, *args):
        base = ["-i", path] if sub == "kaks" else [path]
        a = run_cpp(sub, *base, "-q", "--digits", "17", "-j", "2", *args)
        b = run_cpp(sub, *base, "-q", "--digits", "17", "-j", "2", "--kaks-baseline", *args)
        self.assertEqual(a.returncode, 0, a.stderr)
        self.assertEqual(b.returncode, 0, b.stderr)
        return a.stdout, b.stdout

    def test_same_bytes_on_both_paths(self):
        self.require_avx2_path()
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
        self.require_avx2_path()
        every = ("YN", "NG", "LWL", "MLWL", "LPB", "MLPB", "MYN", "Comeron", "Li")
        fast, base = self.both("dnds", KAKS / "sim_small.axt", *[x for m in every for x in ("-m", m)])
        self.assertEqual(fast, base)

    def test_the_option_is_not_advertised(self):
        self.assertNotIn(b"kaks-baseline", run_cpp("dnds", "--help").stdout)
        self.assertNotIn(b"kaks-baseline", run_cpp("kaks", "--help").stdout)
        self.assertNotIn(b"kaks-isa", run_cpp("dnds", "--help").stdout)
        self.assertNotIn(b"kaks-isa", run_cpp("kaks", "--help").stdout)


if __name__ == "__main__":
    unittest.main()
