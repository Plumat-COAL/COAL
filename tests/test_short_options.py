"""Attached values of one-letter options in the small command lines (cpp/src/shortopts.hpp).

`coal kaks`, `coal dnds`, `coal translate`, `coal codon-dnds` and `coal-cpp-align` read their arguments by hand.
They take `-j4`, `-j=4`, `-mYN`, `-ofile` and chained flags (`-qj4`) as the main command line does (cpp/src/args.cpp,
after Python's argparse), and `--jobs=4` in all five. Each spelling must give the bytes of the spaced one; what is
not an option stays an error; the value of an option is never taken apart.
Skipped if build/coal-cpp (COAL_CPP_BIN) or coal-cpp-align is not built; failed instead with COAL_REQUIRE_CPP.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_align as ta   # noqa: E402  (its lookups of the two binaries)

ROOT = Path(__file__).resolve().parent.parent
AXT = ROOT / "tests" / "kaks" / "sim_small.axt"
CDS = ROOT / "tests" / "fixtures" / "pipeline" / "query_cds.fasta"
GROUPS = "G\tg1\tall\nM\ta\tATGGCTAAAGGTCTGTTTGATCCGGAATGG\nM\tb\tATGGCAAAGGGTTTGTTCGACCCGGAGTGG\nE\n"


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=300, **kw)


class ShortOptionsBase(unittest.TestCase):
    def setUp(self):
        self.cpp = ta.cpp_binary()
        self.align = ta.align_binary()
        if self.cpp is None or self.align is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("coal-cpp / coal-cpp-align not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("coal-cpp / coal-cpp-align not built")
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def same(self, base, *variants, **kw):
        want = run(base, **kw)
        self.assertEqual(want.returncode, 0, want.stderr)
        self.assertTrue(want.stdout)
        for v in variants:
            got = run(v, **kw)
            self.assertEqual((got.returncode, got.stdout), (0, want.stdout), (v, got.stderr))
        return want


class TestAttachedShortOptions(ShortOptionsBase):
    def test_dnds(self):
        c = [self.cpp, "dnds", AXT]
        self.same(c + ["-j", "1", "-q", "-m", "YN", "-m", "Li", "-c", "1"],
                  c + ["-j1", "-q", "-mYN", "-mLi", "-c1"], c + ["-qj1", "-m=YN", "-m=Li", "-c=1"],
                  c + ["-j=1", "-q", "--method=YN", "--method=Li", "--codontable=1"], c + ["-qj=1", "-mYN", "-mLi"],
                  c + ["--jobs=1", "-qmYN", "-mLi"])
        out = self.tmp / "o.tsv"
        r = run(c + ["-qo" + str(out), "-j2"])
        self.assertEqual((r.returncode, r.stdout), (0, ""), r.stderr)
        self.assertEqual(out.read_text(), run(c + ["-q"]).stdout)

    def test_kaks_with_kaks_calculators_letters_in_either_case(self):
        c = [self.cpp, "kaks"]
        self.same(c + ["-i", AXT, "-j", "1", "-q", "-m", "NG", "-c", "1"],
                  c + [f"-i{AXT}", "-j1", "-q", "-mNG", "-c1"], c + [f"-I{AXT}", "-qj1", "-MNG", "-C1"],
                  c + [f"-i={AXT}", "-j=1", "-q", "-m=NG", "-C=1"], c + [f"--in={AXT}", "--jobs=1", "-q", "--method=NG"])

    def test_translate(self):
        c = [self.cpp, "translate", CDS]
        self.same(c + ["-q", "-c", "11"], c + ["-qc11"], c + ["-q", "-c=11"], c + ["-c11", "-q"],
                  c + ["--codontable=11", "-q"], c + ["-q", "-codontable", "11"])
        out = self.tmp / "p.fa"
        r = run(c + [f"-qo{out}"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out.read_text(), run(c + ["-q"]).stdout)

    def test_codon_dnds(self):
        g = self.tmp / "g.tsv"
        g.write_text(GROUPS, encoding="utf-8")
        outs = []
        for k, args in enumerate((["-i", g, "-j", "1", "-m", "YN", "-m", "Comeron", "-c", "1", "-q"],
                                  [f"-i{g}", "-j1", "-mYN", "-mComeron", "-c1", "-q"],
                                  [f"-i={g}", "-qj=1", "-m=YN", "-mC95", "--codontable=1"],
                                  [f"--in={g}", "--jobs=1", "--method=YN", "-qmComeron"])):
            out = self.tmp / f"o{k}.tsv"
            r = run([self.cpp, "codon-dnds", *args, ("-o" + str(out)) if k % 2 else "-o", *([] if k % 2 else [out])])
            self.assertEqual(r.returncode, 0, (args, r.stderr))
            outs.append(out.read_text())
        self.assertEqual(len(outs[0].splitlines()), 3)
        self.assertEqual(set(outs), {outs[0]})

    def test_align(self):
        p = self.tmp / "p.fa"
        p.write_text(">a\nMKVLAAGIVGLW\n>b\nMRVLAGIVGLW\n", encoding="utf-8")
        c = [self.align, p]
        want = self.same(c + ["-j", "1", "-m", "BLOSUM45"], c + ["-j1", "-mBLOSUM45"], c + ["-j=1", "-m=BLOSUM45"],
                         c + ["--jobs=1", "--matrix=BLOSUM45", "--gap-open=10", "--end-gaps=free"])
        out = self.tmp / "a.fa"
        r = run(c + ["-j1", "-mBLOSUM45", f"-o{out}"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out.read_text(), want.stdout)
        self.assertEqual(run([self.align, "--version"]).stdout.split()[0], "coal-cpp-align")

    def test_what_is_not_an_option_is_still_refused(self):
        for args, msg in ((["-qx"], "unrecognized argument: -qx"), (["-q=1"], "unrecognized argument: -q=1"),
                          (["-x1"], "unrecognized argument: -x1"), (["-j"], "-j: expected a number"),
                          (["-jx"], "-j: expected a number"), (["-j="], "-j: expected a number"),
                          (["-j-1"], "-j: expected a number"), (["-mXX"], "XX")):
            for sub in ("dnds", "kaks"):
                r = run([self.cpp, sub, *(["-i"] if sub == "kaks" else []), AXT, *args])
                self.assertEqual(r.returncode, 2, (sub, args, r.stderr))
                self.assertIn(msg, r.stderr, (sub, args))
                self.assertEqual(r.stdout, "")
        r = run([self.align, self.tmp / "none.fa", "-x1"])
        self.assertEqual(r.returncode, 1)
        self.assertIn("unknown option -x1", r.stderr)
        r = run([self.align, self.tmp / "none.fa", "-jx"])
        self.assertEqual(r.returncode, 1)
        self.assertIn("-j: not a number", r.stderr)

    def test_the_value_of_an_option_is_not_taken_apart(self):
        """`-o -q1` writes a file called -q1 (the value of -o), as before; a negative number is no option."""
        r = run([self.cpp, "dnds", AXT, "-q", "-o", "-q1"], cwd=str(self.tmp))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue((self.tmp / "-q1").is_file())
        r = run([self.cpp, "dnds", AXT, "-1"])
        self.assertEqual(r.returncode, 2)
        self.assertIn("unrecognized argument: -1", r.stderr)

    def test_python_command_passes_them_on(self):
        """`coal kaks -j1` / `coal dnds -j1` (the reported case) through the Python command."""
        env = {**os.environ, "COAL_CPP": str(self.cpp), "PYTHONPATH": str(ROOT / "src")}
        want = run([self.cpp, "dnds", AXT, "-q", "-j", "1"]).stdout
        r = run([sys.executable, "-m", "coal", "dnds", AXT, "-q", "-j1"], env=env)
        self.assertEqual((r.returncode, r.stdout), (0, want), r.stderr)
        r = run([sys.executable, "-m", "coal", "kaks", "-i", AXT, "-q", "-j1"], env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout, run([self.cpp, "kaks", "-i", AXT, "-q", "-j", "1"]).stdout)


if __name__ == "__main__":
    unittest.main()
