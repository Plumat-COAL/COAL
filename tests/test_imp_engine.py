"""Improvement list, engine and command line (IMP-E-*, overnight/integration/CHANGES.md).

Run:  python -m unittest tests.test_imp_engine   (C++ checks need COAL_CPP_BIN or ./build/coal-cpp; else skipped)
"""
from __future__ import annotations

import io
import os
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


if __name__ == "__main__":
    unittest.main()
