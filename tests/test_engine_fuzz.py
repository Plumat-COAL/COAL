"""tests/tools/engine_fuzz.py (the nightly coal-cpp vs Python generator) keeps working: a short run passes, a
difference makes it exit 1 and keeps the case, a sanitizer report in coal-cpp's stderr is caught."""
from __future__ import annotations

import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tests" / "tools" / "engine_fuzz.py"
sys.path.insert(0, str(TOOL.parent))
import engine_fuzz  # noqa: E402

CPP = engine_fuzz.cpp_binary()


def run_tool(*args, env=None, timeout=600):
    return subprocess.run([sys.executable, str(TOOL), *args], cwd=ROOT, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=timeout,
                          env=None if env is None else {**os.environ, **env})


class TestGenerator(unittest.TestCase):
    """The generator itself: no engine needed."""

    def test_cases_are_reproducible_and_varied(self):
        import random
        a = [engine_fuzz.make_case(random.Random(5)) for _ in range(3)]
        b = [engine_fuzz.make_case(random.Random(5)) for _ in range(3)]
        self.assertEqual(a, b)
        files, opts, _ = a[0]
        self.assertTrue(any(f.startswith("pep") for f in files))
        self.assertIn("nuc.fasta", files)
        self.assertIn("-output", opts)
        rng = random.Random(11)
        formats = {engine_fuzz.make_case(rng)[1][1] for _ in range(60)}
        self.assertEqual(formats, set(engine_fuzz.FORMATS))

    def test_sanitizer_pattern(self):
        for s in ("==12==ERROR: AddressSanitizer: heap-use-after-free", "x.cpp:3:5: runtime error: signed integer",
                  "==1==ERROR: LeakSanitizer: detected memory leaks"):
            self.assertTrue(engine_fuzz.SANITIZER.search(s), s)
        self.assertFalse(engine_fuzz.SANITIZER.search("coal: error: seq0: protein and DNA are inconsistent"))

    def test_no_binary_exits_2(self):
        r = run_tool("--n", "1", "--cpp", str(ROOT / "no-such-coal-cpp"))
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("coal-cpp not found", r.stderr)


@unittest.skipIf(CPP is None, "coal-cpp not built (build/coal-cpp or COAL_CPP_BIN)")
class TestEngineFuzzTool(unittest.TestCase):
    def test_short_run_passes(self):
        r = run_tool("--n", "20", "--seed", "20261022")
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        self.assertIn("seed 20261022: 20 cases, 0 differences, 0 sanitizer reports", r.stdout)

    def test_sanitizer_mode_runs_cpp_only(self):
        r = run_tool("--n", "5", "--seed", "3", "--sanitizer")
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        self.assertIn("sanitizer run (coal-cpp only), seed 3: 5 cases", r.stdout)

    @unittest.skipIf(os.name == "nt", "uses a POSIX shell wrapper around coal-cpp")
    def test_difference_and_sanitizer_report_fail_and_keep_the_case(self):
        with tempfile.TemporaryDirectory() as t:
            # a coal-cpp whose stderr carries a UBSan line: a difference from Python and a sanitizer report
            wrapper = Path(t) / "coal-cpp"
            wrapper.write_text(f'#!/bin/sh\necho "core.cpp:1:1: runtime error: fake" >&2\nexec "{CPP}" "$@"\n')
            wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
            keep = Path(t) / "keep"
            r = run_tool("--n", "2", "--seed", "1", "--cpp", str(wrapper), "--keep", str(keep))
            self.assertEqual(r.returncode, 1, r.stdout[-2000:] + r.stderr[-2000:])
            self.assertIn("2 differences, 2 sanitizer reports", r.stdout)
            self.assertIn("SANITIZER #0", r.stdout)
            self.assertIn("DIFF #0", r.stdout)
            self.assertTrue((keep / "c0" / "args.txt").is_file())
            self.assertTrue((keep / "san0" / "args.txt").is_file())
            args = (keep / "c0" / "args.txt").read_text().split()
            self.assertTrue((keep / "c0" / args[0]).is_file(), args)  # the protein file named in args.txt
            # --sanitizer alone: the report, not the comparison
            r = run_tool("--n", "1", "--seed", "1", "--cpp", str(wrapper), "--sanitizer")
            self.assertEqual(r.returncode, 1)
            self.assertIn("0 differences, 1 sanitizer reports", r.stdout)


if __name__ == "__main__":
    unittest.main()
