"""coal-cpp command line == `python -m coal` on the reviewed branch (CPP-4).

The reference is the reviewed CLI on origin/keanu_overnight/integration (its src/coal package is
extracted with `git ls-tree` / `git show` into a temp dir, not merged; or COAL_REF_SRC = a src dir).

  * Live comparison (Python 3.10/3.11 only, whose argparse messages coal-cpp reproduces): every
    command line in CASES is run through both; same exit code, same last stderr line on errors, same
    stdout on successful runs. `--help` must equal the reference parser's help after the documented
    C++ changes (HELP_DRIVER: -output choices, -j/-strict help, the C++-only options).
  * Pinned expectations (any Python): the error messages captured from the reference with Python 3.10.
  * C++-only behaviour: -output clustal/paml/codon accepted, --strict "not implemented yet", the C++
    extras match exact names only (so `--j` still means --jobs as in Python).

`python tests/test_cpp_cli.py --write-help` regenerates cpp/src/help_text.hpp from the reference.
Skipped if build/coal-cpp (COAL_CPP_BIN) is not built.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF_BRANCH = "origin/keanu_overnight/integration"
TIMEOUT = 60
LIVE_OK = sys.version_info[:2] in ((3, 10), (3, 11))
P, N = "pep.fasta", "nuc.fasta"


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        # on Windows only .exe files: a shared folder may also hold the Linux build (CPP-8)
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def reference_src(tmp: Path) -> Path | None:
    """A directory containing the reviewed `coal` package (or None if git/the branch is unavailable)."""
    env = os.environ.get("COAL_REF_SRC")
    if env:
        return Path(env)
    # CHANGED(MRG-6): was the integration branch read with git show (before the merge); now this branch's own
    # engine, which is the merged, reviewed reference. See overnight/integration/CHANGES.md
    if (ROOT / "src" / "coal" / "cli.py").is_file():
        return ROOT / "src"
    try:
        r = subprocess.run(["git", "ls-tree", "-r", "--name-only", REF_BRANCH, "src/coal"], cwd=ROOT,
                           capture_output=True, text=True, timeout=TIMEOUT)
        if r.returncode != 0 or not r.stdout.strip():
            return None
        for name in r.stdout.split():
            if not name.endswith(".py"):
                continue
            f = subprocess.run(["git", "show", f"{REF_BRANCH}:{name}"], cwd=ROOT, capture_output=True,
                               timeout=TIMEOUT)
            if f.returncode != 0:
                return None
            dst = tmp / "ref" / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(f.stdout)
    except (OSError, subprocess.SubprocessError):
        return None
    return tmp / "ref" / "src"


# The reference parser with the documented C++ changes, printed with --help formatting.
HELP_DRIVER = r'''
import argparse, sys
import coal.cli as cli

class _Got(Exception):
    pass

def _grab(self, *a, **k):
    raise _Got(self)

argparse.ArgumentParser.parse_args = _grab
try:
    cli.main([])
except _Got as g:
    ap = g.args[0]
for act in ap._actions:
    if act.dest == "output_format":
        act.choices = ["fasta", "clustal", "paml", "codon"]
        act.help = "output format (fasta, clustal, paml or codon)"
    elif act.dest == "jobs":
        act.help = "worker threads (default 0 = all hardware threads)"
    elif act.dest == "strict":
        act.help = "NOT IMPLEMENTED YET in coal-cpp (exits with an error). In the Python engine: " + act.help
g = ap.add_argument_group("C++ engine extras (exact option names only)")
g.add_argument("-keep-going", "--keep-going", action="store_true",
               help="skip sequences that cannot be converted; exit code 2 if any were skipped")
g.add_argument("--max-warnings", type=int, metavar="N",
               help="print at most N per-codon warning lines (default 0 = no limit)")
g.add_argument("--progress", action="store_true",
               help="percent done and time left on stderr (inputs over 10,000 sequences)")
g.add_argument("--resume", action="store_true",
               help="continue an interrupted -o run from its checkpoint (same inputs/options); with -nogap, "
                    "-nomismatch, clustal or codon output only the first pass resumes and the final pass is redone")
g.add_argument("--report", metavar="FILE.json",
               help="write a JSON run report (schema coal.report/1), also when the run fails")
g.add_argument("--json", action="store_true", help="print the run report to stdout (requires -o)")
sys.stdout.write(ap.format_help())
'''


def ref_env(src: Path) -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(src)
    env["COLUMNS"] = "80"
    env.pop("COAL_IMPL", None)
    return env


def expected_help(src: Path) -> str:
    r = subprocess.run([sys.executable, "-c", HELP_DRIVER], env=ref_env(src), capture_output=True, text=True,
                       timeout=TIMEOUT)
    if r.returncode != 0:
        raise AssertionError(f"help driver failed: {r.stderr}")
    return r.stdout


def last_line(text: str) -> str:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    return lines[-1] if lines else ""


# Command lines (run in a temp dir holding pep.fasta and nuc.fasta from tests/fixtures/clean).
ACCEPTED = [
    [P, N], [P, N, "-codontable", "2"], [P, N, "--cod", "2"], [P, N, "-cod", "2"], [P, N, "-codontable=2"],
    [P, N, "--codontable=11"], [P, N, "-codontable", " 1_1 "], [P, N, "-codontable", "+011"],
    [P, N, "-nog"], [P, N, "--nogap"], [P, N, "-nom"], [P, N, "--nomismatch", "-nogap"],
    [P, N, "-j4"], [P, N, "--jobs=2"], [P, N, "-j=3"], [P, N, "--jobs", "1_0"], [P, N, "--j", "1"],
    [P, N, "-q"], [P, N, "-qq"], [P, N, "-qj2"], [P, N, "--qu"],
    ["--", P, N], [P, "--", N], ["-q", "--", P, N], [P, "-o", "x.fa", N], ["-o", "out.fa", P, N],
    [P, N, "-o", "x.fa", "-o", "y.fa"], [P, N, "-output", "fasta"], [P, N, "--output-format=fasta"],
    [P, N, "--output-f", "fasta"], [P, "--", N, "--"],
]
HELP_CASES = [["-h"], ["--help"], ["--h"], [P, "--he"], ["-qh"], ["--help", "-j", "x"], ["-h", "--foo"]]
VERSION_CASES = [["--version"], ["--v"], ["--vers", P]]
ERRORS = [
    [], [P], ["--foo"], ["-q"], [P, "-o", "x.fa"],
    [P, N, "-codontable", "7"], ["-codontable", "x", P, N], [P, N, "-codontable", "2.0"],
    [P, N, "-codontable"], [P, N, "-codontable", "--foo"], [P, N, "--codontable=99999999999999999999"],
    [P, N, "-codontable="], [P, N, "-codontable", "0x2"], [P, N, "-codontable", "-1"],
    [P, N, "-j", "x"], [P, N, "-jx"], [P, N, "-j"], [P, N, "-j", "1__0"], [P, N, "-j", "it's"],
    [P, N, "-j", "a\\b"], [P, N, "-j", 'say "hi"'], [P, N, "-o"], [P, N, "-o", "-q"], [P, N, "--out", "--"],
    [P, N, "-qx"], [P, N, "-qqx"], [P, N, "-q=1"], [P, N, "--quiet=1"], [P, N, "--nogap=1"],
    [P, N, "-nogap=12"], [P, N, "-strict=1"], [P, N, "--strict=yes"], [P, N, "-hx"],
    [P, N, "-out", "x"], [P, N, "-n"], [P, N, "--no"], [P, N, "--o", "x"], [P, N, "-no"],
    [P, N, "--foo"], [P, N, "-x"], [P, N, "-cod=2"], [P, N, "-o", "x.fa", N], [P, N, "-q", N, "-q", N],
    [P, N, "-j", "x", "--help"], [P, N, "-n", "--help"], ["--foo", "-n"],
    [P, N, "-output", "xml"], [P, N, "-output"], [P, N, "--output-format=XML"], [P, N, "-outp", "fasta"],
]
# The C++ engine accepts four formats, so its -output choices list differs (documented).
CHOICES_DIFF = "argument -output/--output-format: invalid choice:"

# Pinned: last stderr line of the reference (Python 3.10) for each ERRORS command line; exit code 2.
PINNED = {
    '': 'coal: error: the following arguments are required: pep_aln, nuc_fasta',
    "'pep.fasta'": 'coal: error: the following arguments are required: nuc_fasta',
    "'--foo'": 'coal: error: the following arguments are required: pep_aln, nuc_fasta',
    "'-q'": 'coal: error: the following arguments are required: pep_aln, nuc_fasta',
    "'pep.fasta' '-o' 'x.fa'": 'coal: error: the following arguments are required: nuc_fasta',
    "'pep.fasta' 'nuc.fasta' '-codontable' '7'": 'coal: error: argument -codontable/--codontable: invalid choice: 7 (choose from 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)',
    "'-codontable' 'x' 'pep.fasta' 'nuc.fasta'": "coal: error: argument -codontable/--codontable: invalid int value: 'x'",
    "'pep.fasta' 'nuc.fasta' '-codontable' '2.0'": "coal: error: argument -codontable/--codontable: invalid int value: '2.0'",
    "'pep.fasta' 'nuc.fasta' '-codontable'": 'coal: error: argument -codontable/--codontable: expected one argument',
    "'pep.fasta' 'nuc.fasta' '-codontable' '--foo'": 'coal: error: argument -codontable/--codontable: expected one argument',
    "'pep.fasta' 'nuc.fasta' '--codontable=99999999999999999999'": 'coal: error: argument -codontable/--codontable: invalid choice: 99999999999999999999 (choose from 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)',
    "'pep.fasta' 'nuc.fasta' '-codontable='": "coal: error: argument -codontable/--codontable: invalid int value: ''",
    "'pep.fasta' 'nuc.fasta' '-codontable' '0x2'": "coal: error: argument -codontable/--codontable: invalid int value: '0x2'",
    "'pep.fasta' 'nuc.fasta' '-codontable' '-1'": 'coal: error: argument -codontable/--codontable: invalid choice: -1 (choose from 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)',
    "'pep.fasta' 'nuc.fasta' '-j' 'x'": "coal: error: argument -j/--jobs: invalid int value: 'x'",
    "'pep.fasta' 'nuc.fasta' '-jx'": "coal: error: argument -j/--jobs: invalid int value: 'x'",
    "'pep.fasta' 'nuc.fasta' '-j'": 'coal: error: argument -j/--jobs: expected one argument',
    "'pep.fasta' 'nuc.fasta' '-j' '1__0'": "coal: error: argument -j/--jobs: invalid int value: '1__0'",
    '\'pep.fasta\' \'nuc.fasta\' \'-j\' "it\'s"': 'coal: error: argument -j/--jobs: invalid int value: "it\'s"',
    "'pep.fasta' 'nuc.fasta' '-j' 'a\\\\b'": "coal: error: argument -j/--jobs: invalid int value: 'a\\\\b'",
    '\'pep.fasta\' \'nuc.fasta\' \'-j\' \'say "hi"\'': 'coal: error: argument -j/--jobs: invalid int value: \'say "hi"\'',
    "'pep.fasta' 'nuc.fasta' '-o'": 'coal: error: argument -o/--out: expected one argument',
    "'pep.fasta' 'nuc.fasta' '-o' '-q'": 'coal: error: argument -o/--out: expected one argument',
    "'pep.fasta' 'nuc.fasta' '--out' '--'": 'coal: error: argument -o/--out: expected one argument',
    "'pep.fasta' 'nuc.fasta' '-qx'": "coal: error: argument -q/--quiet: ignored explicit argument 'x'",
    "'pep.fasta' 'nuc.fasta' '-qqx'": "coal: error: argument -q/--quiet: ignored explicit argument 'x'",
    "'pep.fasta' 'nuc.fasta' '-q=1'": "coal: error: argument -q/--quiet: ignored explicit argument '1'",
    "'pep.fasta' 'nuc.fasta' '--quiet=1'": "coal: error: argument -q/--quiet: ignored explicit argument '1'",
    "'pep.fasta' 'nuc.fasta' '--nogap=1'": "coal: error: argument -nogap/--nogap: ignored explicit argument '1'",
    "'pep.fasta' 'nuc.fasta' '-nogap=12'": "coal: error: argument -nogap/--nogap: ignored explicit argument '12'",
    "'pep.fasta' 'nuc.fasta' '-strict=1'": "coal: error: argument -strict/--strict: ignored explicit argument '1'",
    "'pep.fasta' 'nuc.fasta' '--strict=yes'": "coal: error: argument -strict/--strict: ignored explicit argument 'yes'",
    "'pep.fasta' 'nuc.fasta' '-hx'": "coal: error: argument -h/--help: ignored explicit argument 'x'",
    "'pep.fasta' 'nuc.fasta' '-out' 'x'": 'coal: error: ambiguous option: -out could match -output, -o',
    "'pep.fasta' 'nuc.fasta' '-n'": 'coal: error: ambiguous option: -n could match -nogap, -nomismatch',
    "'pep.fasta' 'nuc.fasta' '--no'": 'coal: error: ambiguous option: --no could match --nogap, --nomismatch',
    "'pep.fasta' 'nuc.fasta' '--o' 'x'": 'coal: error: ambiguous option: --o could match --output-format, --out',
    "'pep.fasta' 'nuc.fasta' '-no'": 'coal: error: ambiguous option: -no could match -nogap, -nomismatch',
    "'pep.fasta' 'nuc.fasta' '--foo'": 'coal: error: unrecognized arguments: --foo',
    "'pep.fasta' 'nuc.fasta' '-x'": 'coal: error: unrecognized arguments: -x',
    "'pep.fasta' 'nuc.fasta' '-cod=2'": 'coal: error: unrecognized arguments: -cod=2',
    "'pep.fasta' 'nuc.fasta' '-o' 'x.fa' 'nuc.fasta'": 'coal: error: unrecognized arguments: nuc.fasta',
    "'pep.fasta' 'nuc.fasta' '-q' 'nuc.fasta' '-q' 'nuc.fasta'": 'coal: error: unrecognized arguments: nuc.fasta nuc.fasta',
    "'pep.fasta' 'nuc.fasta' '-j' 'x' '--help'": "coal: error: argument -j/--jobs: invalid int value: 'x'",
    "'pep.fasta' 'nuc.fasta' '-n' '--help'": 'coal: error: ambiguous option: -n could match -nogap, -nomismatch',
    "'--foo' '-n'": 'coal: error: ambiguous option: -n could match -nogap, -nomismatch',
    "'pep.fasta' 'nuc.fasta' '-output' 'xml'": "coal: error: argument -output/--output-format: invalid choice: 'xml' (choose from 'fasta')",
    "'pep.fasta' 'nuc.fasta' '-output'": 'coal: error: argument -output/--output-format: expected one argument',
    "'pep.fasta' 'nuc.fasta' '--output-format=XML'": "coal: error: argument -output/--output-format: invalid choice: 'XML' (choose from 'fasta')",
    "'pep.fasta' 'nuc.fasta' '-outp' 'fasta'": 'coal: error: ambiguous option: -outp could match -output, -o',
}


def run(cmd: list[str], cwd: Path, env=None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=TIMEOUT)


class _Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bin = cpp_binary()
        if cls.bin is None:
            raise unittest.SkipTest("coal-cpp not built (build/coal-cpp or COAL_CPP_BIN)")
        cls._tmp = tempfile.TemporaryDirectory()
        cls.tmp = Path(cls._tmp.name)
        cls.work = cls.tmp / "work"
        cls.work.mkdir()
        for f in (P, N):
            shutil.copy(ROOT / "tests" / "fixtures" / "clean" / f, cls.work / f)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def cpp(self, args):
        return run([str(self.bin), *args], self.work)


class LiveComparison(_Base):
    @classmethod
    def setUpClass(cls):
        if not LIVE_OK:
            raise unittest.SkipTest("live comparison needs Python 3.10/3.11 (argparse messages differ later)")
        super().setUpClass()
        cls.src = reference_src(cls.tmp)
        if cls.src is None or not (cls.src / "coal" / "cli.py").is_file():
            cls._tmp.cleanup()
            raise unittest.SkipTest(f"reference CLI not available ({REF_BRANCH} or COAL_REF_SRC)")

    def ref(self, args):
        return run([sys.executable, "-m", "coal", *args], self.work, env=ref_env(self.src))

    def test_help_text(self):
        want = expected_help(self.src)
        got = self.cpp(["--help"])
        self.assertEqual(got.returncode, 0)
        self.assertEqual(got.stdout, want)

    def test_errors(self):
        for args in ERRORS:
            with self.subTest(args=args):
                r, c = self.ref(args), self.cpp(args)
                self.assertEqual(r.returncode, 2, r.stderr)
                self.assertEqual(c.returncode, r.returncode, c.stderr)
                rl, cl = last_line(r.stderr), last_line(c.stderr)
                if rl.startswith("coal: error: " + CHOICES_DIFF):
                    rl, cl = rl.split(" (choose from")[0], cl.split(" (choose from")[0]
                self.assertEqual(cl, rl)
                self.assertTrue(c.stderr.startswith("usage: coal "), c.stderr)

    def test_accepted(self):
        for args in ACCEPTED:
            with self.subTest(args=args):
                r, c = self.ref(args), self.cpp(args)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(c.returncode, 0, c.stderr)
                self.assertEqual(c.stdout, r.stdout)

    def test_help_and_version(self):
        for args in HELP_CASES:
            with self.subTest(args=args):
                r, c = self.ref(args), self.cpp(args)
                self.assertEqual((c.returncode, r.returncode), (0, 0))
                self.assertTrue(c.stdout.startswith("usage: coal "), c.stdout[:80])
        for args in VERSION_CASES:
            with self.subTest(args=args):
                r, c = self.ref(args), self.cpp(args)
                self.assertEqual((c.returncode, r.returncode), (0, 0))
                self.assertTrue(c.stdout.startswith("coal "), c.stdout)

    def test_pinned_table_is_current(self):
        # the pinned messages below must still be what the reference prints
        for args in ERRORS:
            with self.subTest(args=args):
                self.assertEqual(PINNED.get(key(args)), last_line(self.ref(args).stderr))


def key(args) -> str:
    return " ".join(repr(a) for a in args)


class Pinned(_Base):
    def test_errors(self):
        self.assertEqual(len(PINNED), len(ERRORS))
        for args in ERRORS:
            with self.subTest(args=args):
                c = self.cpp(args)
                self.assertEqual(c.returncode, 2, c.stderr)
                want, got = PINNED[key(args)], last_line(c.stderr)
                if want.startswith("coal: error: " + CHOICES_DIFF):
                    want, got = want.split(" (choose from")[0], got.split(" (choose from")[0]
                self.assertEqual(got, want)

    def test_accepted(self):
        for args in ACCEPTED:
            with self.subTest(args=args):
                self.assertEqual(self.cpp(args).returncode, 0)


class CppOnly(_Base):
    def test_more_formats_accepted(self):
        for fmt in ("clustal", "paml", "codon"):
            with self.subTest(fmt=fmt):
                c = self.cpp([P, N, "-output", fmt])
                self.assertEqual(c.returncode, 0, c.stderr)
                self.assertTrue(c.stdout)
        c = self.cpp([P, N, "-output", "xml"])
        self.assertEqual(last_line(c.stderr), "coal: error: argument -output/--output-format: invalid choice: "
                                              "'xml' (choose from 'fasta', 'clustal', 'paml', 'codon')")

    def test_strict_not_implemented(self):
        for flag in ("-strict", "--strict", "--str"):
            with self.subTest(flag=flag):
                c = self.cpp([P, N, flag])
                self.assertEqual(c.returncode, 2)
                self.assertEqual(c.stdout, "")
                self.assertIn("not implemented yet", last_line(c.stderr))

    def test_extras_exact_names_only(self):
        c = self.cpp([P, N, "--keep"])
        self.assertEqual(last_line(c.stderr), "coal: error: unrecognized arguments: --keep")
        self.assertEqual(self.cpp([P, N, "--j", "1"]).returncode, 0)  # not ambiguous with --json
        self.assertEqual(self.cpp([P, N, "-keep-going", "--max-warnings=3"]).returncode, 0)
        c = self.cpp([P, N, "--max-warnings", "x"])
        self.assertEqual(last_line(c.stderr), "coal: error: argument --max-warnings: invalid int value: 'x'")
        c = self.cpp([P, N, "--json"])
        self.assertEqual(c.returncode, 2)
        self.assertIn("--json requires -o", last_line(c.stderr))

    def test_version(self):
        c = self.cpp(["--version"])
        self.assertEqual(c.returncode, 0)
        self.assertRegex(c.stdout, r"^coal \d+\.\d+\.\d+")


def _main():
    if "--write-help" in sys.argv or "--print-pinned" in sys.argv:
        with tempfile.TemporaryDirectory() as t:
            src = reference_src(Path(t))
            if src is None:
                sys.exit("reference CLI not available")
            if "--write-help" in sys.argv:
                text = expected_help(src)
                hdr = ROOT / "cpp" / "src" / "help_text.hpp"
                hdr.write_text(
                    "// Generated by `python tests/test_cpp_cli.py --write-help` from the reviewed Python CLI\n"
                    f"// ({REF_BRANCH}) plus the documented C++ changes. Do not edit by hand.\n"
                    "#pragma once\n\n"
                    f'static const char* const COAL_HELP = R"COALHELP({text})COALHELP";\n',
                    encoding="utf-8", newline="\n")
                print(f"wrote {hdr}")
            else:
                work = Path(t) / "work"
                work.mkdir()
                for f in (P, N):
                    shutil.copy(ROOT / "tests" / "fixtures" / "clean" / f, work / f)
                for args in ERRORS:
                    r = run([sys.executable, "-m", "coal", *args], work, env=ref_env(src))
                    print(f"    {key(args)!r}: {last_line(r.stderr)!r},")
        return
    unittest.main()


if __name__ == "__main__":
    _main()
