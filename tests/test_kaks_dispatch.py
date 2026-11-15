"""YN / MYN: the AVX2 copy of the transition-probability matrix gives the same bits as the baseline copy.

On x86 with GCC or Clang, cpp/src/kaks.cpp compiles pmat_codon (the 61 x 61 eigendecomposition of each round of
YN and MYN, most of their time) twice, for the baseline ISA (SSE2) and for AVX2, and uses the AVX2 copy when the CPU
has it (a run-time check; no build or package needs an -m flag). Both copies do the same floating-point operations
in the same order, so the outputs must be byte-identical. `--kaks-baseline` (not in the help; testing only) makes a
run use the baseline copy, and `--kaks-isa` (testing only too) prints which copy a run would use. This test first
checks that the AVX2 copy is the one used (where the CPU has AVX2 and the build has the copy), then compares the two
at 17 significant digits: both output formats, both YN flavours and weightings, several genetic codes. Where there
is no AVX2 copy to compare the comparison is skipped with that reason, so a pass always means both paths ran: an ARM
CPU or a 32-bit x86 build (the copy is compiled for x86-64 only), an MSVC build (no target attribute), a build with
-DCOAL_KAKS_BASELINE (the baseline only, at build time), a CPU without AVX2. Which compiler built coal-cpp, and
whether with COAL_KAKS_BASELINE, is read from its CMake build tree (CMakeFiles/<version>/CMakeCXXCompiler.cmake and
CMakeCache.txt, in the binary's folder or the one above it for build/Release/coal-cpp.exe), else from the binary
itself (GCC's "GCC: (" and Clang's "clang version" strings; the Rich header of an MSVC-linked .exe), not from the
operating system: a MinGW / Rtools GCC build on Windows must have the copy, and one without it fails. A
COAL_KAKS_BASELINE binary away from its build tree cannot be told from a GCC / Clang build that lost the copy: there
three of the four tests fail (test_which_copy_runs, test_same_bytes_on_both_paths, test_every_method_on_both_paths;
test_the_option_is_not_advertised passes). Needs build/coal-cpp (COAL_CPP_BIN).
"""
from __future__ import annotations

import os
import platform
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import CPP, needs_cpp, run_cpp  # noqa: E402
from test_kaks_yn_shortcut import low_divergence_pairs, multi_difference  # noqa: E402

KAKS = ROOT / "tests" / "kaks"
X86 = ("x86_64", "amd64", "i386", "i686", "x86")
X86_64 = ("x86_64", "amd64")   # the AVX2 copy is compiled for x86-64 only (__x86_64__ with SSE2 maths)
CONFIGS = ("Release", "Debug", "RelWithDebInfo", "MinSizeRel")   # a multi-config generator's output folders


def cpu_has_avx2() -> bool | None:
    """True / False where the operating system says, None where this cannot tell."""
    if platform.machine().lower() not in X86:
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
        if os.name == "nt":
            import ctypes
            # PF_AVX2_INSTRUCTIONS_AVAILABLE (40): Windows 10 and later
            return bool(ctypes.windll.kernel32.IsProcessorFeaturePresent(40))
    except (OSError, AttributeError):
        return None
    return None


def build_info(binary) -> tuple[str | None, bool]:
    """(compiler, COAL_KAKS_BASELINE) of the build a coal-cpp comes from: the compiler as CMake names it ("GNU",
    "Clang", "AppleClang", "MSVC", ...; None when nothing says), and whether the build defined COAL_KAKS_BASELINE.
    From the CMake build tree the binary sits in (or, for a multi-config generator's Release/ or Debug/, the one
    above it; only a tree whose CMakeCache.txt names itself), else from marks the compilers leave in the binary.
    clang-cl is "Clang" in both (CMake's id; "clang version" in the binary, looked for before MSVC's Rich header):
    it has the target attribute, so its builds have the copy."""
    binary = Path(binary).resolve()
    trees = [binary.parent] + ([binary.parent.parent] if binary.parent.name in CONFIGS else [])
    for tree in trees:
        cache = tree / "CMakeCache.txt"
        if not cache.is_file():
            continue
        text = cache.read_text(encoding="utf-8", errors="replace")
        home = re.search(r"(?m)^CMAKE_CACHEFILE_DIR:INTERNAL=(.*)$", text)
        if home and Path(home.group(1).strip()).resolve() != tree.resolve():
            continue   # a CMakeCache.txt copied there, not this binary's build tree
        flags = " ".join(m.group(1) for m in re.finditer(r"(?m)^CMAKE_CXX_FLAGS[\w]*:STRING=(.*)$", text))
        baseline = "COAL_KAKS_BASELINE" in flags
        for f in sorted(tree.glob("CMakeFiles/*/CMakeCXXCompiler.cmake")):
            m = re.search(r'(?m)^set\(CMAKE_CXX_COMPILER_ID "([^"]*)"\)', f.read_text(encoding="utf-8", errors="replace"))
            if m and m.group(1):
                return m.group(1), baseline
        return None, baseline
    try:
        data = binary.read_bytes()
    except OSError:
        return None, False
    if b"clang version" in data or b"Apple clang" in data:
        return "Clang", False
    if b"GCC: (" in data:
        return "GNU", False
    if data[:2] == b"MZ" and b"Rich" in data[:1024]:   # the Rich header: written by Microsoft's linker only
        return "MSVC", False
    return None, False


def binary_machine(binary) -> str | None:
    """"x86_64", "x86" (32-bit), "arm64" or None: what the executable itself is built for (ELF, PE; a Mach-O,
    which may hold several, gives None), not what the operating system runs on."""
    try:
        with open(binary, "rb") as fh:
            head = fh.read(4096)
    except OSError:
        return None
    if head[:4] == b"\x7fELF" and len(head) > 20:
        machine = int.from_bytes(head[18:20], "little" if head[5] == 1 else "big")
        return {62: "x86_64", 3: "x86", 183: "arm64"}.get(machine)
    if head[:2] == b"MZ" and len(head) > 0x40:
        pe = int.from_bytes(head[0x3C:0x40], "little")
        if head[pe:pe + 4] == b"PE\0\0":
            machine = int.from_bytes(head[pe + 4:pe + 6], "little")
            return {0x8664: "x86_64", 0x14C: "x86", 0xAA64: "arm64"}.get(machine)
    return None


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
        if platform.machine().lower() not in X86:
            self.skipTest(f"not an x86 CPU ({platform.machine()}): no AVX2 copy is compiled, only the baseline runs")
        if platform.machine().lower() not in X86_64 or binary_machine(CPP) == "x86":
            self.skipTest("a 32-bit x86 build: the AVX2 copy is compiled for x86-64 only, only the baseline runs")
        compiler, baseline = build_info(CPP)
        if compiler == "MSVC":
            self.skipTest("this build has no AVX2 copy (MSVC has no target attribute; cpp/src/kaks.cpp COAL_KAKS_AVX2)")
        if baseline:
            self.skipTest("this build has no AVX2 copy (built with -DCOAL_KAKS_BASELINE: the baseline only)")
        if avx2 is False:
            self.skipTest("this CPU has no AVX2: only the baseline copy can run")
        if avx2 is None:
            self.skipTest("cannot tell whether this CPU has AVX2, and the build runs the baseline copy")
        self.fail(f"the CPU has AVX2 but coal-cpp runs the baseline copy: the AVX2 copy is missing from this "
                  f"{compiler or 'GCC / Clang'} build on x86 (cpp/src/kaks.cpp COAL_KAKS_AVX2)")

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
