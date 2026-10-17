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
# CHANGED(MRG-6): was None (skipped) without COAL_REF_SRC; after the merge this branch's own engine is the reference
REF_SRC = Path(os.environ["COAL_REF_SRC"]).resolve() if os.environ.get("COAL_REF_SRC") else ROOT / "src"
TIMEOUT = 60


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def run_cpp(args, cwd, env=None):
    r = subprocess.run([str(cpp_binary()), *args], cwd=cwd, capture_output=True, timeout=TIMEOUT,
                       env=None if env is None else dict(os.environ, **env))
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


SHIM_CLI = """\"\"\"coal-cpp in place of the reference's `coal` command (built by tests/test_cpp_vs_reference.py).\"\"\"
import os
import subprocess
import sys

CPP = {cpp!r}


def _on_sigterm(signum, frame):
    raise SystemExit(128 + signum)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    sys.stdout.flush()
    sys.stderr.flush()
    return subprocess.call([CPP, *argv], env=dict(os.environ, COAL_CPP_ONLY_LINES="0"))
"""
SHIM_MAIN = """import os
import sys

CPP = {cpp!r}
os.environ["COAL_CPP_ONLY_LINES"] = "0"  # leave out the stderr lines only coal-cpp prints
if os.name == "nt":
    import subprocess
    sys.exit(subprocess.call([CPP, *sys.argv[1:]]))
os.execv(CPP, [CPP, *sys.argv[1:]])
"""
SHIM_INIT = """import os as _os
__path__.append({ref!r})  # every module but cli / __main__ comes from the reference package
with open(_os.path.join({ref!r}, "__init__.py")) as _f:
    exec(compile(_f.read(), "coal/__init__.py", "exec"))
"""
REF_TEST_MODULES = ("test_decisions", "test_review_fixes", "test_fuzz_findings", "test_final_fixes")


def build_shim(tmp: Path) -> Path:
    """A copy of the reference's test layout whose `coal` command (python -m coal, and coal.cli.main used by the
    `pal2nal` drop-in) runs coal-cpp; every other module is the reference's. Returns the tests folder."""
    ref_root = REF_SRC.parent
    src = tmp / "src"
    (src / "coal").mkdir(parents=True)
    cpp = str(cpp_binary())
    (src / "coal" / "__init__.py").write_text(SHIM_INIT.format(ref=str(REF_SRC / "coal")))
    (src / "coal" / "__main__.py").write_text(SHIM_MAIN.format(cpp=cpp))
    (src / "coal" / "cli.py").write_text(SHIM_CLI.format(cpp=cpp))
    for name in ("coal_compat", "coal_app"):
        if (REF_SRC / name).exists():
            os.symlink(REF_SRC / name, src / name)
    (tmp / "tests").mkdir()
    for name in ("fixtures", "golden", "cases.json"):
        os.symlink(ref_root / "tests" / name, tmp / "tests" / name)
    for name in ("scripts", "benchmarks", "examples", "README.md", "pyproject.toml"):
        if (ref_root / name).exists():
            os.symlink(ref_root / name, tmp / name)
    for m in REF_TEST_MODULES:
        (tmp / "tests" / f"{m}.py").write_bytes((ref_root / "tests" / f"{m}.py").read_bytes())
    return tmp / "tests"


def run_ref_tests(modules=REF_TEST_MODULES, timeout=170):
    import tempfile
    with tempfile.TemporaryDirectory(prefix="coal_shim_") as t:
        tests = build_shim(Path(t))
        env = dict(os.environ, PYTHONPATH=str(Path(t) / "src"))
        r = subprocess.run([sys.executable, "-m", "unittest", "-v", *modules], cwd=tests, env=env,
                           capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout + r.stderr


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


# Reference tests that cannot apply to a separate binary: they replace sys.stdout inside the Python process.
NOT_APPLICABLE = ("TestDec3ClosedPipe.test_posix_broken_pipe_is_141_and_silent_in_both_commands",
                  "TestDec3ClosedPipe.test_windows_einval_is_141_too")


@unittest.skipIf(_skip_reason(), _skip_reason() or "")
class TestReferenceTestsWithCpp(unittest.TestCase):
    """The reference's own engine tests (test_decisions, test_review_fixes, test_fuzz_findings), with its `coal`
    command (python -m coal, and the `pal2nal` drop-in's engine call) running coal-cpp."""

    def test_reference_tests(self):
        rc, out = run_ref_tests()
        bad = [ln for ln in out.splitlines() if ln.endswith(("... FAIL", "... ERROR"))
               and not any(n.split(".")[1] in ln and n.split(".")[0] in ln for n in NOT_APPLICABLE)]
        self.assertTrue(" tests in " in out, out[-2000:])
        self.assertEqual(bad, [], out[-4000:])


@unittest.skipIf(_skip_reason(), _skip_reason() or "")
class TestCppVsReferenceLarge(unittest.TestCase):
    """CPP-20 / CPP-24: 1,200 sequences (PAML header widths; FASTA and marked CLUSTAL input), every format
    and every clustal/codon writer path, coal-cpp vs the reference Python engine."""

    def test_many_sequences_every_format(self):
        import random
        import tempfile
        rng = random.Random(11)
        code = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
                "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
                "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            pep, nuc = [], []
            for k in range(1200):
                prot = "M" + "".join(rng.choice(sorted(code)) for _ in range(29))
                pep.append(f">s{k}\n{prot}\n")
                nuc.append(f">s{k}\n{''.join(code[a] for a in prot)}\n")
            (d / "p.fa").write_text("".join(pep))
            (d / "n.fa").write_text("".join(nuc))
            # CPP-24: the same rows as a CLUSTAL file with '#' block marks (mark lines under clustal / codon blocks)
            rows = [(f"s{k}", p.split("\n")[1]) for k, p in enumerate(pep)]
            aln = ["CLUSTAL W\n\n"]
            for b in range(0, 30, 12):
                aln += [f"{i:<10}{r[b:b + 12]}\n" for i, r in rows]
                aln.append(" " * 10 + "".join(rng.choice(" #") for _ in range(min(12, 30 - b))) + "\n\n")
            (d / "p.aln").write_text("".join(aln))
            env = dict(os.environ, PYTHONPATH=str(REF_SRC))
            # writer variants: default (chunked, one chunk), many chunks, the few-blocks-of-every-row writer
            variants = {"default": {}, "chunks": {"COAL_BLOCK_CHUNK": "40000"}, "rows": {"COAL_BLOCK_CHUNK": "1"}}
            for pep_file in ("p.fa", "p.aln"):
                for fmt in ("paml", "fasta", "clustal", "codon"):
                    ref = subprocess.run([sys.executable, "-m", "coal", pep_file, "n.fa", "-output", fmt], cwd=d,
                                         env=env, capture_output=True, timeout=TIMEOUT)
                    self.assertEqual(ref.returncode, 0, ref.stderr[-500:])
                    self.assertGreater(len(ref.stdout), 100000)
                    if pep_file == "p.aln" and fmt in ("clustal", "codon"):
                        self.assertIn(b"#", ref.stdout)  # the mark lines are there
                    for name, venv in variants.items():
                        with self.subTest(pep=pep_file, fmt=fmt, writer=name):
                            rc, out, _ = run_cpp([pep_file, "n.fa", "-output", fmt], d, venv)
                            self.assertEqual((rc, out), (ref.returncode, ref.stdout.decode("latin-1")))


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
    elif what == "reftests":
        rc, out = run_ref_tests(tuple(sys.argv[2:]) or REF_TEST_MODULES)
        print(out[-int(os.environ.get("SHOW", "20000")):])
    else:
        unittest.main()
