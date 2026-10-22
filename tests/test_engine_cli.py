"""Engine and command line: logged finish errors, inert test hooks, acceptance sets, genetic-code hints, the same
stderr from both engines, the warning cap, help groups and automatic coal-cpp dispatch.

Run:  python -m unittest tests.test_engine_cli   (C++ checks need COAL_CPP_BIN or ./build/coal-cpp; else skipped)
"""
from __future__ import annotations

import io
import os
import re
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal import report as report_mod  # noqa: E402
from coal.cli import main as coal_main  # noqa: E402

FIX = ROOT / "tests" / "fixtures"


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe")):
        if cand.is_file():
            return cand
    return None


def run_py(args, cwd=None):
    out, err = io.StringIO(), io.StringIO()
    old = os.getcwd()
    if cwd:
        os.chdir(cwd)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = coal_main(list(args))
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(old)
    return rc, out.getvalue(), err.getvalue()


class TestFinishErrorsAreLogged(unittest.TestCase):
    """IMP-E-G: an exception in the report's finish step is printed under COAL_DEBUG; the run goes on."""

    def test_debug_prints_and_run_goes_on(self):
        d = FIX / "codon_mismatch"
        args = [str(d / "pep.fasta"), str(d / "nuc.fasta")]
        rc0, out0, _ = run_py(args)
        with mock.patch.object(report_mod.RunReport, "_finish", side_effect=RuntimeError("boom")):
            with mock.patch.dict(os.environ, {"COAL_DEBUG": "1"}):
                rc, out, err = run_py(args)
            self.assertEqual((rc, out), (rc0, out0))
            self.assertIn("coal: debug: the report's finish step failed", err)
            self.assertIn("RuntimeError: boom", err)
            with mock.patch.dict(os.environ, {"COAL_DEBUG": ""}):
                rc, out, err = run_py(args)
            self.assertEqual((rc, out), (rc0, out0))
            self.assertNotIn("debug", err)


def make_synthetic(d: Path, n: int, codons: int, seed: int = 1) -> tuple[str, str]:
    """benchmarks/gen_data.py's model (standard code); returns (pep, nuc) paths."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from gen_data import generate
    generate(n, codons, d, seed)
    return str(d / "pep.fasta"), str(d / "nuc.fasta")


class TestTestHooksInert(unittest.TestCase):
    """IMP-E-I: COAL_TEST_PAUSE_* do nothing unless COAL_TEST_HOOKS=1 is set too."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coal_hooks_"))
        cls.pep, cls.nuc = make_synthetic(cls.tmp, 300, 100)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _run(self, cmd):
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_BATCH_NT": "20000", "COAL_CHECKPOINT_EVERY": "0",
               "COAL_TEST_PAUSE_AFTER_CHECKPOINT": "1", "COAL_TEST_PAUSE_BEFORE_CHECKPOINT": "1",
               "COAL_ENGINE": "python"}
        env.pop("COAL_TEST_HOOKS", None)
        r = subprocess.run(cmd + [self.pep, self.nuc, "-o", str(self.tmp / "out.fa")], env=env,
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("coal: test:", r.stderr)

    def test_python(self):
        self._run([sys.executable, "-m", "coal"])

    def test_cpp(self):
        b = cpp_binary()
        if b is None:
            self.skipTest("coal-cpp not built")
        self._run([str(b)])


class TestFastAcceptanceSets(unittest.TestCase):
    """IMP-E-H: accepted() and _note_ambiguous's quick test give what the brute-force versions gave."""

    def test_accepted_equals_brute_force(self):
        from itertools import product
        from coal import codon
        from coal.tables import SUPPORTED
        for table in SUPPORTED:
            for aa in "ACDEFGHIKLMNPQRSTVWY*":
                for mode in (codon.NORMAL, codon.START, codon.ANCHOR_START):
                    allowed = codon._allowed(table, aa, mode)
                    brute = {"".join(t) for t in product(codon.IUPAC, repeat=3)
                             if all("".join(c) in allowed for c in product(*(codon.IUPAC[b] for b in t)))}
                    self.assertEqual(codon.accepted(table, aa, mode), brute, (table, aa, mode))

    def test_quick_residue_test_equals_counting(self):
        from coal import core
        self.assertEqual(core.AMBIGUOUS_MAX_RESIDUES, 3)  # the regex counts to 4
        rng = random.Random(7)
        for _ in range(20000):
            row = "".join(rng.choice("--------0123456789ACMW*") for _ in range(rng.randint(0, 12)))
            old = len(row) - row.count("-") - sum(row.count(d) for d in "0123456789") > 3
            self.assertEqual(bool(core._MORE_THAN_3_RES.search(row)), old, row)


class TestCheapHint(unittest.TestCase):
    """IMP-E-D: the hint is skipped under -q (unless --report/--json) and tests at most HINT_SAMPLE sequences."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="coal_hint_"))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.pep, self.nuc = make_synthetic(self.tmp, 120, 60)

    def _calls(self, extra):
        import coal.core
        with mock.patch.object(coal.core, "thread_one", wraps=coal.core.thread_one) as m:
            rc, out, err = run_py([self.pep, self.nuc, "-codontable", "2", *extra])
        return rc, err, m.call_count

    def test_quiet_skips_the_hint(self):
        rc, err, n = self._calls(["-q"])
        self.assertEqual((rc, err, n), (0, "", 0))
        rc, err, n = self._calls(["-q", "-o", str(self.tmp / "o.fa"), "--report", str(self.tmp / "r.json")])
        self.assertGreater(n, 0)
        import json
        h = json.loads((self.tmp / "r.json").read_text())["hints"]
        self.assertEqual([x["suggested_codontable"] for x in h if x["code"] == "genetic_code"], [1])

    def test_sample_and_early_stop(self):
        rc, err, n = self._calls([])
        self.assertIn(f"COAL tested {report_mod.HINT_SAMPLE} of the 120 sequences", err)
        self.assertIn("Re-run with -codontable 1.", err)
        # every sampled gene under code 2 and under each code that ties for the best (1, 4, 6, 10, 11, 15, 16 here),
        # the other 9 codes only until they can no longer win: 431 calls, not 17 x 50 = 850 (the early notice of
        # IMP-E-B tests the same 50 first; the end-of-run hint reuses its answers)
        self.assertLess(n, 9 * report_mod.HINT_SAMPLE)


class TestSuggestionFirst(unittest.TestCase):
    """IMP-E-B: one `coal: notice:` with the genetic-code suggestion before the per-codon warnings, naming the code the
    end-of-run summary names."""

    def _check(self, args, code, cwd=None):
        rc, out, err = run_py(args, cwd)
        lines = err.splitlines()
        first = next(i for i, x in enumerate(lines) if "pepAlnPos" in x)
        notes = [i for i, x in enumerate(lines) if x.startswith("coal: notice: wrong genetic code?")]
        self.assertEqual(len(notes), 1, err)
        self.assertLess(notes[0], first)
        self.assertTrue(lines[notes[0]].endswith(f"Re-run with -codontable {code}."), lines[notes[0]])
        self.assertTrue(lines[-1].endswith(f"Re-run with -codontable {code}."), lines[-1])
        return lines

    def test_mito_example_with_code_1(self):
        lines = self._check(["protein.aln.fasta", "dna.fasta"], 2, cwd=ROOT / "examples" / "mito")
        self.assertIn("all 24 sequences match their protein; under code 1, none do.", lines[0])

    def test_synthetic_with_code_2(self):
        tmp = Path(tempfile.mkdtemp(prefix="coal_first_"))
        self.addCleanup(shutil.rmtree, tmp, True)
        pep, nuc = make_synthetic(tmp, 80, 60)
        lines = self._check([pep, nuc, "-codontable", "2"], 1)
        self.assertIn("50 of the first 50 sequences match their protein; under code 2, none do.", lines[0])

    def test_no_notice_when_quiet_nomismatch_or_right_code(self):
        d = ROOT / "examples" / "mito"
        for extra in (["-q"], ["-nomismatch"], ["-codontable", "2"]):
            rc, out, err = run_py(["protein.aln.fasta", "dna.fasta", *extra], d)
            self.assertNotIn("wrong genetic code?", err, extra)


@unittest.skipIf(cpp_binary() is None, "coal-cpp not built (set COAL_CPP_BIN)")
class TestStderrSameBothEngines(unittest.TestCase):
    """IMP-E-E: on wrong-code runs both engines print the same stderr byte for byte (notice, warnings, summary)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coal_e_"))
        sys.path.insert(0, str(ROOT / "tests" / "realdata"))
        import realdata_lib as L
        rows = L.clean_rows(L.load(ROOT / "tests" / "realdata" / "coxI_mammals" / "COXI.nex.gz", "nexus"))
        pep, nuc, _ = L.derive(rows, 2)
        (cls.tmp / "cox").mkdir()
        (cls.tmp / "cox" / "pep.fasta").write_text(pep)
        (cls.tmp / "cox" / "nuc.fasta").write_text(nuc)
        cls.syn = make_synthetic(cls.tmp / "syn", 300, 120)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _both(self, args, cwd):
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_ENGINE": "python"}
        py = subprocess.run([sys.executable, "-m", "coal", *args], cwd=cwd, env=env, capture_output=True, timeout=120)
        cp = subprocess.run([str(cpp_binary()), *args], cwd=cwd, env=env, capture_output=True, timeout=120)
        return py, cp

    def test_wrong_code_runs(self):
        cox = self.tmp / "cox"
        runs = [(["pep.fasta", "nuc.fasta"], cox), (["pep.fasta", "nuc.fasta", "-codontable", "11"], cox),
                (["protein.aln.fasta", "dna.fasta"], ROOT / "examples" / "mito"),
                (["protein.aln.fasta", "dna.fasta", "-codontable", "2"], ROOT / "examples" / "one_broken_gene"),
                ([*self.syn, "-codontable", "2"], self.tmp)]
        for args, cwd in runs:
            for extra in ([], ["--max-warnings", "0"], ["--max-warnings", "3"], ["-q"], ["--keep-going"]):
                with self.subTest(args=args, extra=extra):
                    py, cp = self._both(args + extra, cwd)
                    self.assertEqual(py.returncode, cp.returncode)
                    self.assertEqual(py.stdout, cp.stdout)
                    self.assertEqual(py.stderr.decode().replace("\r\n", "\n"), cp.stderr.decode().replace("\r\n", "\n"))
                    if not extra or extra == ["--keep-going"]:
                        self.assertIn(b"Re-run with -codontable", py.stderr)


class TestWarningCap(unittest.TestCase):
    """IMP-E-C: `coal` caps per-codon warnings at 20 by default; the `pal2nal` drop-in prints them all."""

    def test_default_cap_and_drop_in(self):
        d = ROOT / "examples" / "mito"
        rc, out, err = run_py(["protein.aln.fasta", "dna.fasta"], d)
        rc0, out0, err0 = run_py(["protein.aln.fasta", "dna.fasta", "--max-warnings", "0"], d)
        n = err0.count("pepAlnPos")
        self.assertGreater(n, 20)
        self.assertEqual((rc, out), (rc0, out0))
        self.assertEqual(err.count("pepAlnPos"), 20)
        self.assertIn(f"coal: warning: ... and {n - 20:,} more warnings (see them all with --max-warnings 0 or "
                      f"--report FILE.json)", err)
        from coal_compat.cli import main as pal2nal_main
        o, e = io.StringIO(), io.StringIO()
        old = os.getcwd()
        os.chdir(d)
        try:
            with redirect_stdout(o), redirect_stderr(e):
                pal2nal_main(["protein.aln.fasta", "dna.fasta", "-output", "fasta"])
        finally:
            os.chdir(old)
        self.assertEqual(e.getvalue().count("pepAlnPos"), n)
        self.assertNotIn("more warnings", e.getvalue())


class TestHelpGroups(unittest.TestCase):
    """IMP-E-F: --help in four groups with short help; coal-cpp prints the same help (one documented line, -j) and
    accepts every option spelling the Python parser has."""

    def test_groups_and_length(self):
        from coal.cli import build_parser
        h = build_parser().format_help()
        for g in ("Core:", "PAL2NAL options:", "Large runs:", "Reporting:"):
            self.assertIn("\n" + g + "\n", h)
        self.assertLess(len(h.splitlines()), 60)
        # the app reads the formats here (coal_app.core.detect_features). Python <= 3.12 prints
        # "-output {...}, --output-format {...}", so the check uses the app's own pattern.
        m = re.search(r"-output(?:,\s*--[\w-]+)?\s+\{([\w,]+)\}", h)
        self.assertIsNotNone(m, h)
        self.assertEqual(m.group(1), "fasta,clustal,paml,codon")

    @unittest.skipIf(cpp_binary() is None, "coal-cpp not built")
    # coal-cpp's help mirrors Python 3.13's argparse layout ("-j, --jobs N"); 3.12 and older print
    # "-j N, --jobs N", so the line-for-line comparison only holds on 3.13+ (CI's cpp job). Spellings are checked below.
    @unittest.skipIf(sys.version_info < (3, 13), "help layout differs before Python 3.13")
    def test_cpp_help_and_spellings(self):
        from coal.cli import build_parser
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COLUMNS": "80", "COAL_ENGINE": "python"}
        py = subprocess.run([sys.executable, "-m", "coal", "--help"], env=env, capture_output=True, text=True,
                            timeout=60).stdout.splitlines()
        cp = subprocess.run([str(cpp_binary()), "--help"], env=env, capture_output=True, text=True,
                            timeout=60).stdout.splitlines()
        i, j = next(k for k, x in enumerate(py) if "-j, --jobs" in x), next(k for k, x in enumerate(cp) if "-j, --jobs" in x)
        self.assertEqual(py[:i] + py[i + 2:], cp[:j] + cp[j + 1:])  # -j: "worker threads" in C++, one line
        self.assertIn("worker threads", cp[j])
        values = {"-codontable": "2", "--codontable": "2", "-output": "fasta", "--output-format": "fasta",
                  "-o": "x.fa", "--out": "x.fa", "-j": "1", "--jobs": "1", "--max-warnings": "3",
                  "--report": "r.json"}
        d = FIX / "clean"
        tmp = Path(tempfile.mkdtemp(prefix="coal_sp_"))
        self.addCleanup(shutil.rmtree, tmp, True)
        for act in build_parser()._actions:
            for opt in act.option_strings:
                if opt in ("-h", "--help", "--version"):
                    continue
                args = [str(d / "pep.fasta"), str(d / "nuc.fasta"), opt] + ([values[opt]] if opt in values else [])
                if opt in ("--json", "--resume"):
                    args += ["-o", "y.fa"]
                with self.subTest(opt=opt):
                    p = subprocess.run([sys.executable, "-m", "coal", *args], cwd=tmp, env=env, capture_output=True,
                                       timeout=60)
                    c = subprocess.run([str(cpp_binary()), *args], cwd=tmp, env=env, capture_output=True, timeout=60)
                    self.assertEqual(p.returncode, 0, p.stderr)
                    self.assertEqual(c.returncode, 0, c.stderr)


@unittest.skipIf(os.name == "nt", "uses a POSIX shell wrapper to see which engine ran")
@unittest.skipIf(cpp_binary() is None, "coal-cpp not built (set COAL_CPP_BIN)")
class TestAutomaticCppDispatch(unittest.TestCase):
    """IMP-E-A: `coal` runs coal-cpp for large inputs, with the same stdout, stderr, report and exit code."""

    @classmethod
    def setUpClass(cls):
        from coal import dispatch
        cls.tmp = Path(tempfile.mkdtemp(prefix="coal_disp_"))
        cls.big = make_synthetic(cls.tmp / "big", 600, 500)  # ~1.17 MB > CPP_MIN_BYTES
        cls.small = make_synthetic(cls.tmp / "small", 20, 100)
        assert sum(os.path.getsize(p) for p in cls.big) >= dispatch.CPP_MIN_BYTES
        assert sum(os.path.getsize(p) for p in cls.small) < dispatch.CPP_MIN_BYTES
        cls.mark = cls.tmp / "cpp_ran"
        cls.wrapper = cls.tmp / "coal-cpp"
        cls.wrapper.write_text(f'#!/bin/sh\ncase "$1" in --version) ;; *) touch "{cls.mark}" ;; esac\n'
                               f'exec "{cpp_binary()}" "$@"\n')
        cls.wrapper.chmod(0o755)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_coal(self, args, engine=None, cpp=None, module="coal"):
        if self.mark.exists():
            self.mark.unlink()
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_CPP": str(cpp or self.wrapper)}
        env.pop("COAL_ENGINE", None)
        if engine:
            env["COAL_ENGINE"] = engine
        work = self.tmp / "work"
        if work.exists():
            shutil.rmtree(work)
        work.mkdir()
        r = subprocess.run([sys.executable, "-m", module, *args], cwd=work, env=env, capture_output=True, timeout=300)
        rep = (work / "r.json").read_text() if (work / "r.json").exists() else None
        out = (work / "o.fa").read_bytes() if (work / "o.fa").exists() else None
        import json
        if rep:
            rep = json.loads(rep)
            rep.pop("runtime_seconds")
        stdout = r.stdout
        if "--json" in args:  # the report on stdout: the same but for its run time
            stdout = json.loads(stdout)
            stdout.pop("runtime_seconds")
        return r.returncode, stdout, r.stderr, out, rep, self.mark.exists()

    def test_large_input_uses_cpp_with_identical_results(self):
        for extra in ([], ["-codontable", "2"], ["-codontable", "2", "-o", "o.fa", "--report", "r.json"],
                      ["-output", "clustal", "-nogap", "-q"], ["-codontable", "2", "--max-warnings", "0"],
                      ["-o", "o.fa", "--json"]):
            with self.subTest(extra=extra):
                auto = self.run_coal([*self.big, *extra])
                py = self.run_coal([*self.big, *extra], engine="python")
                self.assertTrue(auto[5], "coal-cpp did not run")
                self.assertFalse(py[5])
                self.assertEqual(auto[:5], py[:5])

    def test_small_input_stays_in_python(self):
        auto = self.run_coal([*self.small, "-codontable", "2"])
        self.assertFalse(auto[5])
        self.assertIn(b"Re-run with -codontable 1.", auto[2])

    def test_forced_modes(self):
        r = self.run_coal([*self.small, "-codontable", "2"], engine="cpp")
        self.assertTrue(r[5])
        self.assertEqual(r[:5], self.run_coal([*self.small, "-codontable", "2"], engine="python")[:5])
        self.assertFalse(self.run_coal(list(self.big), engine="python")[5])
        rc, out, err, *_ = self.run_coal(list(self.small), engine="cpp", cpp=self.tmp / "missing")
        self.assertEqual((rc, out), (1, b""))
        self.assertIn(b"coal: error: COAL_ENGINE=cpp: coal-cpp was not found", err)
        rc, out, err, *_ = self.run_coal(list(self.small), engine="fast")
        self.assertEqual(rc, 1)
        self.assertIn(b"COAL_ENGINE='fast'", err)

    def test_cpp_that_cannot_start_falls_back_with_one_notice(self):
        broken = self.tmp / "broken-coal-cpp"
        broken.write_text("not a program\n")
        broken.chmod(0o644)
        r = self.run_coal(list(self.big), cpp=broken)
        py = self.run_coal(list(self.big), engine="python")
        self.assertEqual((r[0], r[1]), (py[0], py[1]))
        extra = r[2].decode().replace(py[2].decode(), "", 1)
        self.assertEqual(extra.count("\n"), 1)
        self.assertTrue(extra.startswith("coal: notice: could not start"), extra)
        self.assertTrue(extra.rstrip().endswith("; using the Python engine"), extra)
        stale = self.tmp / "stale-coal-cpp"
        stale.write_text('#!/bin/sh\necho "coal 0.0.1 (C++)"\n')
        stale.chmod(0o755)
        r = self.run_coal(list(self.big), cpp=stale)
        self.assertEqual((r[0], r[1]), (py[0], py[1]))
        self.assertIn(b"is coal 0.0.1 (C++), not coal ", r[2])
        self.assertNotIn(b"notice", self.run_coal([*self.big, "-q"], cpp=stale)[2])

    def test_pal2nal_drop_in_stays_in_python(self):
        r = self.run_coal([*self.big, "-output", "fasta"], module="coal_compat")
        self.assertEqual(r[0], 0)
        self.assertFalse(r[5])

    def test_app_gets_dispatch_for_free(self):
        sys.path.insert(1, str(ROOT / "apps" / "coal-app" / "src"))  # the app's source, as report_vs_app.py does
        try:
            from coal_app import core
        except ImportError as e:
            self.skipTest(f"coal_app not importable: {e}")
        if self.mark.exists():
            self.mark.unlink()
        out = self.tmp / "app_out"
        out.mkdir(exist_ok=True)
        env = {**os.environ, "COAL_CPP": str(self.wrapper)}
        env.pop("COAL_ENGINE", None)
        with mock.patch.dict(os.environ, env, clear=True):
            r = core.run_coal(*self.big, codontable=2, out_dir=str(out))
        self.assertTrue(self.mark.exists(), "the app's `python -m coal` did not run coal-cpp")
        self.assertEqual([h["suggested_codontable"] for h in r.report["hints"] if h["code"] == "genetic_code"], [1])

    def test_binary_lookup_order(self):
        from coal import dispatch
        with mock.patch.dict(os.environ, {"COAL_CPP": str(self.wrapper)}):
            self.assertEqual(dispatch.find_cpp(), str(self.wrapper))
        with mock.patch.dict(os.environ, {"COAL_CPP": str(self.tmp / "missing")}):
            self.assertIsNone(dispatch.find_cpp())
        fake_py = self.tmp / "venv" / "bin" / "python"
        fake_py.parent.mkdir(parents=True, exist_ok=True)
        internal = self.tmp / "venv" / "bin" / "_internal"
        internal.mkdir(exist_ok=True)
        shutil.copy(self.wrapper, internal / "coal-cpp")
        env = {k: v for k, v in os.environ.items() if k != "COAL_CPP"}
        with mock.patch.dict(os.environ, env, clear=True), mock.patch.object(sys, "executable", str(fake_py)):
            self.assertEqual(dispatch.find_cpp(), str(internal / "coal-cpp"))  # a frozen app's _internal
            shutil.copy(self.wrapper, fake_py.parent / "coal-cpp")
            self.assertEqual(dispatch.find_cpp(), str(fake_py.parent / "coal-cpp"))  # next to Python first
            os.unlink(fake_py.parent / "coal-cpp")
            os.unlink(internal / "coal-cpp")
            with mock.patch.object(dispatch.sysconfig, "get_path", lambda name: str(self.tmp / "none")), \
                    mock.patch.dict(os.environ, {"PATH": str(self.tmp)}):
                self.assertEqual(dispatch.find_cpp(), str(self.tmp / "coal-cpp"))  # then PATH


if __name__ == "__main__":
    unittest.main()
