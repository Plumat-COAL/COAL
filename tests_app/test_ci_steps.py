"""R2-D4 / R2-D11: run CI step scripts for real instead of grepping them."""
import os
import re
import shutil
import subprocess
import tempfile
import textwrap
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ci():
    with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
        return fh.read()


def step_script(job: str, step_name: str) -> str:
    """The `run: |` script of one step (plain text parsing: CI's environment has no YAML library)."""
    ci = _ci()
    m = re.search(rf"^  {re.escape(job)}:\n(.*?)(?=^  [\w-]+:\n|\Z)", ci, re.S | re.M)
    assert m, job
    s = re.search(rf"^(\s*)- name: {re.escape(step_name)}\n(.*?)(?=^\1- name:|\Z)", m.group(1), re.S | re.M)
    assert s, step_name
    r = re.search(r"^(\s*)run: \|\n((?:\1  .*\n|\s*\n)+)", s.group(2), re.M)
    assert r, f"no run: | block in {step_name}"
    return textwrap.dedent(r.group(2))


BASH = shutil.which("bash")


def _run_cpp_step(tc, step, files=(), python_out=None, python_rc=0):
    """Run one cpp-job step script in a temp dir with the GitHub files set; `python` is a stub printing python_out.
    Returns (CompletedProcess, {"output": ..., "env": ..., "summary": ...})."""
    d = tempfile.mkdtemp(prefix="coal_ci_step_")
    tc.addCleanup(shutil.rmtree, d, True)
    for rel, text in files:
        p = os.path.join(d, *rel.split("/"))
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.chmod(p, 0o755)
    stub = os.path.join(d, "_stub")
    os.makedirs(stub)
    with open(os.path.join(d, "_py_out.txt"), "w", encoding="utf-8") as fh:
        fh.write(python_out or "")
    with open(os.path.join(stub, "python"), "w", encoding="utf-8") as fh:
        fh.write(f'#!/bin/sh\ncat "{d}/_py_out.txt"\nexit {python_rc}\n')
    os.chmod(os.path.join(stub, "python"), 0o755)
    gh = {k: os.path.join(d, f"_gh_{k}") for k in ("output", "env", "summary")}
    env = dict(os.environ, PATH=stub + os.pathsep + os.environ.get("PATH", ""), GITHUB_OUTPUT=gh["output"],
               GITHUB_ENV=gh["env"], GITHUB_STEP_SUMMARY=gh["summary"])
    r = subprocess.run([BASH, "--noprofile", "--norc", "-eo", "pipefail", "-c", step_script("cpp", step)], cwd=d,
                       env=env, capture_output=True, text=True, timeout=60)
    out = {}
    for k, p in gh.items():
        out[k] = open(p, encoding="utf-8").read() if os.path.exists(p) else ""
    return r, out


CPP_TESTS = ["tests/test_differential.py", "tests/test_cpp_vs_python.py", "tests/test_cpp_cli.py",
             "tests/test_cpp_pool.py", "tests/test_cpp_readers.py", "tests/test_cpp_tables.py",
             "tests/test_cpp_vs_reference.py"]


# CHANGED(CI-CPP1): was TestCppRequireStep for cpp/build/coal(.exe) and test_cpp_vs_python.py only, with no test of
# the detect step (gated on cpp/CMakeLists.txt, which never exists). Now the detect, require, golden and C++ test
# steps of the job as it builds the root CMakeLists.txt. See overnight/integration/CHANGES.md
@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash")
class TestCppDetectStep(unittest.TestCase):  # CI-CPP1
    STEP = "Check for the root CMakeLists.txt with coal-cpp"

    def test_repo_root_cmakelists_is_detected(self):
        with open(os.path.join(ROOT, "CMakeLists.txt"), encoding="utf-8") as fh:
            real = fh.read()
        r, out = _run_cpp_step(self, self.STEP, [("CMakeLists.txt", real)])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("present=true", out["output"])

    def test_missing_or_other_project_is_not_applicable(self):
        for files in ([], [("cpp/CMakeLists.txt", "add_executable(coal-cpp main.cpp)\n")],
                      [("CMakeLists.txt", "add_executable(coal-cpp-dump dump.cpp)\n")]):
            with self.subTest(files=files):
                r, out = _run_cpp_step(self, self.STEP, files)
                self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                self.assertIn("present=false", out["output"])
                self.assertIn("not applicable", out["summary"])


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash")
class TestCppRequireStep(unittest.TestCase):  # R2-D4, CI-CPP1
    STEP = "Require the binary and its C++ tests"

    def run_step(self, binary=None, tests=CPP_TESTS):
        files = [(t, "") for t in tests] + ([(binary, "#!/bin/sh\n")] if binary else [])
        return _run_cpp_step(self, self.STEP, files)

    def test_makefile_binary_passes_and_is_exported(self):
        r, out = self.run_step("build/coal-cpp")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("COAL_CPP_BIN=build/coal-cpp\n", out["env"])

    def test_visual_studio_binary_passes_and_is_exported(self):
        r, out = self.run_step("build/Release/coal-cpp.exe")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("COAL_CPP_BIN=build/Release/coal-cpp.exe\n", out["env"])

    def test_no_binary_fails(self):
        for binary in (None, "cpp/build/coal"):  # the old location no longer counts
            with self.subTest(binary=binary):
                r, out = self.run_step(binary)
                self.assertEqual(r.returncode, 1)
                self.assertIn("no build/coal-cpp", r.stdout + r.stderr)
                self.assertEqual(out["env"], "")

    def test_missing_cpp_test_module_fails(self):
        for missing in ("tests/test_differential.py", "tests/test_cpp_vs_python.py", "tests/test_cpp_readers.py"):
            with self.subTest(missing=missing):
                r, _ = self.run_step("build/coal-cpp", [t for t in CPP_TESTS if t != missing])
                self.assertEqual(r.returncode, 1)
                self.assertIn(missing, r.stdout + r.stderr)

    def test_listed_modules_exist(self):
        for t in CPP_TESTS:
            self.assertTrue(os.path.isfile(os.path.join(ROOT, *t.split("/"))), t)


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash")
class TestCppGoldenStep(unittest.TestCase):  # CI-CPP1
    STEP = "Golden checks with the C++ engine (must not skip)"
    OK = "cpp: 285/285 golden checks\nPAL2NAL comparison: 39 cases, 285 checks: golden rc = 0 285\n"

    def test_all_passing_passes_and_writes_summary(self):
        r, out = _run_cpp_step(self, self.STEP, python_out=self.OK)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("cpp: 285/285 golden checks", out["summary"])

    def test_skipped_fails(self):
        r, _ = _run_cpp_step(self, self.STEP, python_out="cpp: skipped (coal-cpp not built; ...)\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("did not run against coal-cpp", r.stdout + r.stderr)

    def test_python_engine_count_fails(self):  # COAL_IMPL not honoured
        r, _ = _run_cpp_step(self, self.STEP, python_out="python: 285/285 golden checks\n")
        self.assertEqual(r.returncode, 1)

    def test_failing_checks_fail(self):
        r, _ = _run_cpp_step(self, self.STEP, python_out="clean/fasta: x\ncpp: 284/285 golden checks\n", python_rc=1)
        self.assertNotEqual(r.returncode, 0)

    def test_step_sets_coal_impl_cpp(self):
        self.assertRegex(_ci().split(f"- name: {self.STEP}")[1].split("run: |")[0], r"COAL_IMPL: cpp")


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash")
class TestCppTestsStep(unittest.TestCase):  # CI-CPP1
    STEP = "C++ engine tests (must not skip)"

    def test_runs_every_cpp_module(self):
        script = step_script("cpp", self.STEP)
        for t in CPP_TESTS:
            self.assertIn("tests." + os.path.basename(t)[:-3], script)

    def test_by_design_skips_pass(self):
        out = ("setUpClass (tests.test_cpp_cli.LiveComparison) ... skipped 'live comparison needs Python 3.10/3.11 "
               "(argparse messages differ later)'\ntest_fuzz ... skipped 'set PAL2NAL=/path/to/pal2nal.pl'\n"
               "test_x ... skipped 'POSIX signals'\nRan 40 tests\nOK (skipped=3)\n")
        r, out_ = _run_cpp_step(self, self.STEP, python_out=out)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("Ran 40 tests", out_["summary"])

    def test_engine_not_found_skips_fail(self):
        for msg in ("coal-cpp not built", "coal-cpp not built (set COAL_CPP_BIN or build into ./build)",
                    "coal-cpp-dump not built (cmake --build build)",
                    "reference reader not available (origin/x not fetched; or set COAL_REF_SEQIO)",
                    "set COAL_REF_SRC=<reference worktree>/src"):
            with self.subTest(msg=msg):
                r, _ = _run_cpp_step(self, self.STEP, python_out=f"test_a ... skipped '{msg}'\nRan 1 test\nOK (skipped=1)\n")
                self.assertEqual(r.returncode, 1)
                self.assertIn("were skipped", r.stdout + r.stderr)

    def test_failures_fail(self):
        r, _ = _run_cpp_step(self, self.STEP, python_out="Ran 9 tests\nFAILED (failures=1)\n", python_rc=1)
        self.assertNotEqual(r.returncode, 0)


class TestCppJobShape(unittest.TestCase):  # CI-CPP1 (plain text: no YAML library in CI)
    def setUp(self):
        m = re.search(r"^  cpp:\n(.*?)(?=^  [\w-]+:\n|\Z)", _ci(), re.S | re.M)
        self.job = m.group(1)

    def test_three_os_matrix_with_bash(self):
        self.assertIn("runs-on: ${{ matrix.os }}", self.job)
        self.assertIn("os: [ubuntu-latest, windows-latest, macos-latest]", self.job)
        self.assertRegex(self.job, r"defaults:\n\s+run:\n\s+shell: bash")
        self.assertIn('python-version: "3.13"', self.job)

    def test_builds_root_project_with_warnings_as_errors(self):
        self.assertIn("cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOAL_WARNINGS_AS_ERRORS=ON", self.job)
        self.assertIn("cmake --build build --config Release", self.job)
        code = "\n".join(ln.split(" #")[0] for ln in self.job.splitlines() if not ln.lstrip().startswith("#"))
        self.assertNotIn("cpp/build", code)
        self.assertNotIn("cpp/CMakeLists.txt", code)
        with open(os.path.join(ROOT, "CMakeLists.txt"), encoding="utf-8") as fh:
            self.assertIn("option(COAL_WARNINGS_AS_ERRORS", fh.read())

    def test_name_is_not_stale(self):
        name = re.search(r"^    name: (.*)$", self.job, re.M).group(1).split("  #")[0]
        self.assertIn("coal-cpp", name)
        self.assertNotIn("not applicable", name)


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash (the step runs on ubuntu-latest)")
class TestFuzzStep(unittest.TestCase):  # DEC-4
    STEP = "Fuzz COAL against PAL2NAL (must not skip)"

    def run_step(self, body):
        d = tempfile.mkdtemp(prefix="coal_ci_step_")
        self.addCleanup(shutil.rmtree, d, True)
        os.makedirs(os.path.join(d, "tests"))
        with open(os.path.join(d, "tests", "fuzz_vs_pal2nal.py"), "w", encoding="utf-8") as fh:
            fh.write(body)
        env = dict(os.environ, GITHUB_STEP_SUMMARY=os.path.join(d, "summary.md"))
        script = step_script("pal2nal-fuzz", self.STEP).replace("python ", f"{__import__('sys').executable} ")
        return subprocess.run([BASH, "--noprofile", "--norc", "-eo", "pipefail", "-c", script], cwd=d, env=env,
                              capture_output=True, text=True, timeout=60)

    def test_passing_fuzz_passes_and_writes_summary(self):
        r = self.run_step("import sys\nprint('fuzz vs PAL2NAL v14, -output fasta: 1 cases', file=sys.stderr)\n")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_skipped_fuzz_fails(self):
        r = self.run_step("import sys\nprint('OK (skipped=1)', file=sys.stderr)\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("skipped", r.stdout + r.stderr)

    def test_failing_fuzz_fails(self):
        r = self.run_step("raise SystemExit(1)\n")
        self.assertNotEqual(r.returncode, 0)

    def test_fetch_step_runs_the_fetch_script(self):
        ci = _ci()
        job = ci.split("  pal2nal-fuzz:")[1].split("\n  explainer:")[0]
        self.assertIn("run: python scripts/fetch_pal2nal.py", job)
        self.assertIn("PAL2NAL: PAL2NAL/pal2nal.v14/pal2nal.pl", job)


if __name__ == "__main__":
    unittest.main()
