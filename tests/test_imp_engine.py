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


if __name__ == "__main__":
    unittest.main()
