"""coal-cpp against the REVIEWED reference (CPP-9): the reference's PAL2NAL golden checks and its fuzz against the
real PAL2NAL v14, with coal-cpp in place of the Python engine.

The reference is the reviewed Python engine (origin/keanu_overnight/integration) checked out as a separate folder,
e.g. `git worktree add ../coal-ref origin/keanu_overnight/integration`; COAL_REF_SRC points at its src/ (its tests/
are next to it). Its harness code is loaded as is; only the "run COAL" function is swapped for one that runs
coal-cpp (COAL_CPP_BIN, default build/coal-cpp):
  * golden: the reference's tests/test_differential.py compare_all (every stored PAL2NAL golden, its DELIBERATE
    rules, warnings compared as there);
  * fuzz: the reference's tests/fuzz_vs_pal2nal.py (needs PAL2NAL=/path/to/pal2nal.pl and perl; FUZZ_N, FUZZ_SEED,
    FUZZ_FORMATS as there). The PAL2NAL-style command line is turned into coal-cpp's by the reference's own
    `pal2nal` drop-in translation (coal_compat.cli.translate).
Skipped when COAL_REF_SRC or the binary is missing.

  python tests/test_cpp_vs_reference.py golden      # list every failing golden check
  python tests/test_cpp_vs_reference.py fuzz        # fuzz summary + the first differences
"""
from __future__ import annotations

import contextlib
import importlib.util
import os
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_SRC = Path(os.environ["COAL_REF_SRC"]).resolve() if os.environ.get("COAL_REF_SRC") else None
TIMEOUT = 60


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def run_cpp(args, cwd):
    r = subprocess.run([str(cpp_binary()), *args], cwd=cwd, capture_output=True, timeout=TIMEOUT)
    err = r.stderr.decode("utf-8", "replace").replace("\r\n", "\n")
    err = "".join(ln for ln in err.splitlines(True) if not ln.startswith("coal: summary:"))  # coal-cpp's own lines
    return r.returncode, r.stdout.decode("latin-1"), err


def load_ref(name: str):
    """A module from the reference's tests/ folder, with the reference src first on sys.path."""
    if str(REF_SRC) not in sys.path:
        sys.path.insert(0, str(REF_SRC))
    for m in [m for m in sys.modules if m == "coal" or m.startswith(("coal.", "coal_compat"))]:
        if not str(getattr(sys.modules[m], "__file__", "")).startswith(str(REF_SRC)):
            del sys.modules[m]  # a coal package from another src must not be reused
    path = REF_SRC.parent / "tests" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"ref_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _Recorder:
    """Stands in for the TestCase compare_all expects: records each failing subTest instead of stopping."""

    def __init__(self):
        self._tc = unittest.TestCase()
        self.failures = []

    @contextlib.contextmanager
    def subTest(self, **kw):
        try:
            yield
        except AssertionError as e:
            self.failures.append((kw, str(e).splitlines()[0][:300] if str(e) else "assertion failed"))

    def __getattr__(self, name):
        return getattr(self._tc, name)


def golden_failures():
    td = load_ref("test_differential")

    def run_coal(case, opts):
        spec = td.MANIFEST["cases"][case]
        args = [spec["pep"], *spec.get("nuc", ["nuc.fasta"]), *opts]
        return run_cpp(args, td.FIX / case)

    td.run_coal = run_coal
    rec = _Recorder()
    ncases, n, count = td.compare_all(rec)
    return n, rec.failures, td.summary_line(ncases, n, count)


def fuzz_run(pal2nal, n, seed, formats, report=print):
    fz = load_ref("fuzz_vs_pal2nal")
    from coal_compat.cli import translate, UsageError
    from coal.pipeline import FORMATS as ENGINE_FORMATS

    def run_coal(d, args):
        try:
            cargs = translate(list(args), ENGINE_FORMATS)
        except UsageError as e:
            return 2, "", f"pal2nal: error: {e}\n"
        return run_cpp(cargs, d)

    fz.run_coal = run_coal
    return fz.fuzz(pal2nal, n, seed, report=report, keep_dir=os.environ.get("FUZZ_KEEP"), formats=formats)


def _skip_reason():
    if REF_SRC is None or not (REF_SRC / "coal" / "cli.py").is_file():
        return "set COAL_REF_SRC=<reference worktree>/src"
    if cpp_binary() is None:
        return "coal-cpp not built"
    return None


@unittest.skipIf(_skip_reason(), _skip_reason() or "")
class TestCppVsReferenceGolden(unittest.TestCase):
    def test_golden(self):
        n, failures, summary = golden_failures()
        print("\n" + summary, file=sys.stderr)
        self.assertGreater(n, 200)
        self.assertEqual(failures, [], f"{len(failures)} failing golden checks, first: {failures[:3]}")


@unittest.skipIf(_skip_reason() or not os.environ.get("PAL2NAL"),
                 _skip_reason() or "set PAL2NAL=/path/to/pal2nal.pl (reference scripts/fetch_pal2nal.py)")
class TestCppVsPal2nalFuzz(unittest.TestCase):
    def test_fuzz(self):
        n = int(os.environ.get("FUZZ_N", "1000"))
        seed = int(os.environ.get("FUZZ_SEED", "20261017"))
        formats = tuple(os.environ.get("FUZZ_FORMATS", "clustal,paml,fasta,codon").split(","))
        counts, bad = fuzz_run(os.path.abspath(os.environ["PAL2NAL"]), n, seed, formats,
                               report=lambda s: print(s, file=sys.stderr))
        self.assertEqual(bad, [], f"{len(bad)} unexplained differences, first: {bad[:1]}")


if __name__ == "__main__":
    if _skip_reason():
        sys.exit(_skip_reason())
    what = sys.argv[1] if len(sys.argv) > 1 else "golden"
    if what == "golden":
        n, failures, summary = golden_failures()
        for kw, msg in failures:
            print(f"{kw.get('case')}/{kw.get('variant')}: {msg}")
        print(summary)
        print(f"cpp vs reference goldens: {n - len(failures)}/{n} checks pass")
    elif what == "fuzz":
        n = int(os.environ.get("FUZZ_N", "1000"))
        counts, bad = fuzz_run(os.path.abspath(os.environ["PAL2NAL"]), n,
                               int(os.environ.get("FUZZ_SEED", "20261017")),
                               tuple(os.environ.get("FUZZ_FORMATS", "clustal,paml,fasta,codon").split(",")))
        for b in bad[:int(os.environ.get("FUZZ_SHOW", "8"))]:
            print({k: b[k] for k in ("case", "args", "pal2nal_ok", "coal_rc", "pal2nal_stdout", "coal_stdout",
                                     "coal_stderr")})
        print(f"unexplained differences: {len(bad)}")
    else:
        unittest.main()
