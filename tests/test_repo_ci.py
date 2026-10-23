"""The CI workflow: step scripts run for real (R2-D4 / R2-D11) instead of grepped, required jobs, pinned actions, minimal permissions, and the extras jobs. The nightly workflow (nightly.yml: fuzz vs PAL2NAL with both engines, the engine generator, sanitizers) under the same rules."""
import os
import re
import shutil
import subprocess
import tempfile
import textwrap
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ci(name="ci.yml"):
    with open(os.path.join(ROOT, ".github", "workflows", name), encoding="utf-8") as fh:
        return fh.read()


def step_script(job: str, step_name: str, workflow: str = "ci.yml") -> str:
    """The `run: |` script of one step (plain text parsing: CI's environment has no YAML library)."""
    ci = _ci(workflow)
    m = re.search(rf"^  {re.escape(job)}:\n(.*?)(?=^  [\w-]+:\n|\Z)", ci, re.S | re.M)
    assert m, job
    s = re.search(rf"^(\s*)- name: {re.escape(step_name)}\n(.*?)(?=^\1- name:|\Z)", m.group(1), re.S | re.M)
    assert s, step_name
    r = re.search(r"^(\s*)run: \|\n((?:\1  .*\n|\s*\n)+)", s.group(2), re.M)
    assert r, f"no run: | block in {step_name}"
    return textwrap.dedent(r.group(2))


BASH = shutil.which("bash")


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash (the step runs on ubuntu-latest)")
class TestCppRequireStep(unittest.TestCase):  # R2-D4
    STEP = "Require the binary and its comparison tests"

    def run_step(self, binary=None, with_tests=True):
        d = tempfile.mkdtemp(prefix="coal_ci_step_")
        self.addCleanup(shutil.rmtree, d, True)
        if binary:
            p = os.path.join(d, *binary.split("/"))
            os.makedirs(os.path.dirname(p))
            with open(p, "w", encoding="utf-8") as fh:
                fh.write("#!/bin/sh\n")
            os.chmod(p, 0o755)
        if with_tests:
            os.makedirs(os.path.join(d, "tests"))
            open(os.path.join(d, "tests", "test_cpp_vs_python.py"), "w").close()
        script = step_script("cpp", self.STEP)
        return subprocess.run([BASH, "--noprofile", "--norc", "-eo", "pipefail", "-c", script], cwd=d,
                              capture_output=True, text=True, timeout=60)

    def test_linux_binary_passes(self):
        r = self.run_step("build/coal-cpp")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_windows_binary_passes(self):
        r = self.run_step("build/Release/coal-cpp.exe")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_no_binary_fails(self):
        r = self.run_step(None)
        self.assertEqual(r.returncode, 1)
        self.assertIn("no build/coal-cpp binary", r.stdout + r.stderr)

    def test_missing_comparison_tests_fail(self):
        r = self.run_step("build/coal-cpp", with_tests=False)
        self.assertEqual(r.returncode, 1)
        self.assertIn("test_cpp_vs_python.py", r.stdout + r.stderr)


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


class TestCiConditions(unittest.TestCase):
    def setUp(self):
        import re
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
            self.ci = fh.read()
        self.re = re

    def _job(self, name):
        m = self.re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.ci, self.re.S | self.re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def _step(self, job, name):
        m = self.re.search(rf"^(\s*)- name: {self.re.escape(name)}\n(.*?)(?=^\1- name:|\Z)", self._job(job),
                           self.re.S | self.re.M)
        self.assertIsNotNone(m, name)
        return m.group(2)

    def test_desktop_step_runs_even_after_app_tests_fail(self):
        step = self._step("python", "Desktop wrapper tests")
        cond = self.re.search(r"^\s*if: (.*)$", step, self.re.M)
        self.assertIsNotNone(cond)
        self.assertIn("!cancelled()", cond.group(1))
        self.assertIn("matrix.os != 'macos-latest'", cond.group(1))
        self.assertNotIn("always()", cond.group(1))

    def test_explainer_job_is_required(self):
        """The job and its test step are required; only the best-effort apt step may continue on error (the test
        step fails loudly when Chromium cannot start, so a slow mirror never hides a real failure)."""
        job = self._job("explainer")
        self.assertNotIn("continue-on-error", job.split("steps:")[0])  # not at the job level
        for step in job.split("\n      - ")[1:]:
            if "continue-on-error" in step:
                self.assertIn("Chromium's system libraries", step.splitlines()[0], step)

    # the GUI job and the desktop step are required on every OS, Windows included
    def test_windows_gui_job_is_required(self):
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._job("gui"), self.re.M))

    def test_desktop_step_is_required_on_windows(self):
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._step("python", "Desktop wrapper tests"), self.re.M))


class TestCiHygiene(unittest.TestCase):
    def test_ci_permissions_and_pinned_actions(self):
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
            ci = fh.read()
        self.assertRegex(ci, r"(?m)^permissions:\n  contents: read")
        uses = re.findall(r"uses: (\S+)", ci)
        self.assertTrue(uses)
        for u in uses:
            self.assertRegex(u, r"@[0-9a-f]{40}$", f"unpinned action {u}")
        with open(os.path.join(ROOT, ".github", "workflows", "release.yml"), encoding="utf-8") as fh:
            rel = {r.split("@")[0]: r for r in re.findall(r"uses: (\S+)", fh.read())}
        for u in uses:
            self.assertEqual(u, rel.get(u.split("@")[0], u), "not the SHA release.yml uses")


class TestNightlyWorkflow(unittest.TestCase):
    """nightly.yml: scheduled plus manual with a seed, the same hygiene as ci.yml (pinned SHAs identical to ci.yml's,
    read-only token, a timeout on every job, one run at a time), the three checks, and artifacts on failure."""

    SCHEDULE = "17 3 * * *"

    def setUp(self):
        self.wf = _ci("nightly.yml")
        self.ci = _ci()

    def _job(self, name):
        m = re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.wf, re.S | re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def _jobs(self):
        return re.findall(r"^  ([\w-]+):\n", self.wf.split("\njobs:\n", 1)[1], re.M)

    def test_triggers_seed_permissions_and_concurrency(self):
        self.assertIn(f"cron: '{self.SCHEDULE}'", self.wf)
        self.assertIn("workflow_dispatch:", self.wf)
        self.assertRegex(self.wf, r"(?m)^      seed:\n")  # the optional seed input
        self.assertIn("inputs.seed || github.run_number", self.wf)  # seed: the input, else the run number
        self.assertRegex(self.wf, r"(?m)^run-name: .*inputs\.seed \|\| github\.run_number")  # ... in the title
        self.assertRegex(self.wf, r"(?m)^  FUZZ_SEED: \$\{\{ inputs\.seed \|\| github\.run_number \}\}")
        self.assertRegex(self.wf, r"(?m)^permissions:\n  contents: read")
        self.assertRegex(self.wf, r"(?m)^concurrency:\n  group: nightly-")
        self.assertNotIn("pull_request", self.wf)
        self.assertNotRegex(self.wf, r"(?m)^  push:")

    def test_every_job_has_a_timeout_and_no_write_permission(self):
        jobs = self._jobs()
        self.assertEqual(sorted(jobs), ["engine-fuzz", "pal2nal-fuzz", "sanitizers"])
        for j in jobs:
            job = self._job(j)
            self.assertRegex(job, r"(?m)^    timeout-minutes: \d+", j)
            self.assertNotIn("permissions:", job, j)
            self.assertNotIn("continue-on-error", job, j)

    def test_actions_pinned_to_the_shas_ci_uses(self):
        uses = re.findall(r"uses: (\S+)", self.wf)
        self.assertTrue(uses)
        ci = {u.split("@")[0]: u for u in re.findall(r"uses: (\S+)", self.ci)}
        for u in uses:
            self.assertRegex(u, r"@[0-9a-f]{40}$", f"unpinned action {u}")
            self.assertEqual(u, ci.get(u.split("@")[0]), f"{u}: not the SHA ci.yml uses")

    def test_pal2nal_fuzz_runs_both_engines_with_the_seed_and_keeps_cases(self):
        job = self._job("pal2nal-fuzz")
        self.assertIn("run: python scripts/fetch_pal2nal.py", job)
        self.assertIn("PAL2NAL: PAL2NAL/pal2nal.v14/pal2nal.pl", job)
        self.assertEqual(job.count('FUZZ_N: "3000"'), 2)
        self.assertIn("COAL_FUZZ_ENGINE: cpp", job)
        self.assertIn("FUZZ_KEEP: fuzz-keep/python", job)
        self.assertIn("FUZZ_KEEP: fuzz-keep/cpp", job)
        self.assertIn("if: ${{ !cancelled() }}", job)  # the coal-cpp half runs after a failed Python half
        self.assertNotIn("FUZZ_SEED:", job)  # the workflow-level seed, not a fixed one
        self.assertIn("if: failure()", job)
        self.assertIn("actions/upload-artifact@", job)
        self.assertIn("fuzz-keep/", job)

    def test_engine_generator_job(self):
        job = self._job("engine-fuzz")
        self.assertIn("ubuntu-latest", job)
        self.assertIn("macos-latest", job)
        self.assertIn('tests/tools/engine_fuzz.py --n 1500 --seed "$FUZZ_SEED" --keep engine-fuzz-keep', job)
        self.assertIn("if: failure()", job)
        self.assertIn("engine-fuzz-keep/", job)
        self.assertTrue(os.path.isfile(os.path.join(ROOT, "tests", "tools", "engine_fuzz.py")))
        self.assertTrue(os.path.isfile(os.path.join(ROOT, "tests", "test_engine_fuzz.py")))  # the tool's own test

    def test_sanitizer_job(self):
        job = self._job("sanitizers")
        self.assertIn("-fsanitize=address,undefined -fno-omit-frame-pointer -O1 -g", job)
        self.assertIn('-DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"', job)
        self.assertIn("-DCOAL_WARNINGS_AS_ERRORS=ON", job)
        self.assertIn("ASAN_OPTIONS: detect_leaks=1", job)
        self.assertIn("UBSAN_OPTIONS: print_stacktrace=1:halt_on_error=1", job)
        self.assertIn('COAL_SKIP_MEMORY_TESTS: "1"', job)
        self.assertIn("COAL_IMPL: cpp", job)
        for mod in ("tests.test_differential", "tests.test_cpp_vs_python", "tests.test_cpp_cli", "tests.test_cpp_readers",
                    "tests.test_cpp_tables", "tests.test_cpp_pool", "tests.test_stdin", "tests.test_robustness",
                    "tests.test_cpp_error_precedence", "tests.test_batch"):
            self.assertIn(mod, job, mod)
        self.assertIn('engine_fuzz.py --n 300 --seed "$FUZZ_SEED" --sanitizer', job)
        self.assertEqual(job.count('grep -qE "AddressSanitizer|LeakSanitizer|runtime error"'), 2)

    def test_the_ci_cpp_test_set_is_the_nightly_one(self):
        """ci.yml's C++ job and the sanitizer job run the same core modules (the nightly adds four)."""
        ci_mods = set(re.findall(r"tests\.test_\w+", self.ci.split("  cpp:")[1]))
        ci_mods.discard("tests.test_wheels")  # packaging (the platform wheel installs): no sanitizer business
        san_mods = set(re.findall(r"tests\.test_\w+", self._job("sanitizers")))
        self.assertTrue(ci_mods <= san_mods, ci_mods - san_mods)


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash (the steps run on ubuntu-latest)")
class TestNightlySteps(unittest.TestCase):
    """The nightly step scripts, run for real (as TestFuzzStep) with stand-in scripts."""

    def _run(self, job, step, files, env=None):
        d = tempfile.mkdtemp(prefix="coal_nightly_step_")
        self.addCleanup(shutil.rmtree, d, True)
        for name, body in files.items():
            os.makedirs(os.path.dirname(os.path.join(d, name)), exist_ok=True)
            with open(os.path.join(d, name), "w", encoding="utf-8") as fh:
                fh.write(body)
        e = dict(os.environ, GITHUB_STEP_SUMMARY=os.path.join(d, "summary.md"), FUZZ_SEED="123", MATRIX_OS="test-os",
                 **(env or {}))
        script = re.sub(r"(?m)^(\s*)python ", lambda m: m.group(1) + __import__("sys").executable + " ",
                        step_script(job, step, "nightly.yml"))  # the command only, not tests.test_cpp_vs_python
        return subprocess.run([BASH, "--noprofile", "--norc", "-eo", "pipefail", "-c", script], cwd=d, env=e,
                              capture_output=True, text=True, timeout=60), d

    PY_STEP = "Fuzz the Python engine against PAL2NAL (must not skip)"
    CPP_STEP = "Fuzz coal-cpp against PAL2NAL (must not skip)"
    FUZZ = "tests/fuzz_vs_pal2nal.py"

    def test_python_fuzz_step_passes_and_writes_the_seed_to_the_summary(self):
        r, d = self._run("pal2nal-fuzz", self.PY_STEP,
                         {self.FUZZ: "import sys\nprint('fuzz vs PAL2NAL v14, -output fasta: 1 cases', file=sys.stderr)\n"})
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(d, "summary.md"), encoding="utf-8") as fh:
            self.assertIn("seed 123", fh.read())

    def test_fuzz_steps_fail_on_a_skip_or_a_failure(self):
        for step in (self.PY_STEP, self.CPP_STEP):
            r, _ = self._run("pal2nal-fuzz", step, {self.FUZZ: "import sys\nprint('OK (skipped=1)', file=sys.stderr)\n"})
            self.assertEqual(r.returncode, 1, step)
            self.assertIn("skipped", r.stdout + r.stderr)
            r, _ = self._run("pal2nal-fuzz", step, {self.FUZZ: "raise SystemExit(1)\n"})
            self.assertNotEqual(r.returncode, 0, step)

    def test_cpp_fuzz_step_requires_the_cpp_engine_line(self):
        # Under `-v` unittest has already printed "test_... ... " on that line, so the marker is not at a line start
        # (the first real nightly run failed on exactly that)
        ok = ("import sys\nprint('test_random_inputs_identical (...) ... ', end='', file=sys.stderr)\n"
              "print('fuzz engine cpp, seed 123, 1 cases', file=sys.stderr)\n")
        r, _ = self._run("pal2nal-fuzz", self.CPP_STEP, {self.FUZZ: ok})
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        r, _ = self._run("pal2nal-fuzz", self.CPP_STEP, {self.FUZZ: ok.replace("cpp", "python")})
        self.assertEqual(r.returncode, 1)
        self.assertIn("did not run coal-cpp", r.stdout + r.stderr)

    def test_sanitizer_test_step_fails_on_a_report_even_if_the_tests_pass(self):
        step = "C++ tests under the sanitizers"
        files = {"tests/__init__.py": ""}
        for mod in ("test_differential", "test_cpp_vs_python", "test_cpp_cli", "test_cpp_readers", "test_cpp_tables",
                    "test_cpp_pool", "test_stdin", "test_robustness", "test_cpp_error_precedence", "test_batch"):
            files[f"tests/{mod}.py"] = "import unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n        pass\n"
        r, _ = self._run("sanitizers", step, files)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        leaky = dict(files)
        leaky["tests/test_batch.py"] = ("import sys, unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n"
                                        "        print('==7==ERROR: LeakSanitizer: detected memory leaks', file=sys.stderr)\n")
        r, _ = self._run("sanitizers", step, leaky)
        self.assertEqual(r.returncode, 1)
        self.assertIn("a sanitizer reported", r.stdout + r.stderr)
        unbuilt = dict(files)
        unbuilt["tests/test_batch.py"] = ("import unittest\nclass T(unittest.TestCase):\n"
                                          "    @unittest.skip('coal-cpp not built')\n    def test_x(self):\n        pass\n")
        r, _ = self._run("sanitizers", step, unbuilt)
        self.assertEqual(r.returncode, 1)
        self.assertIn("not built", r.stdout + r.stderr)

    def test_generator_steps_fail_when_the_tool_fails(self):
        tool = "tests/tools/engine_fuzz.py"
        passing = {tool: "print('engine fuzz, coal-cpp vs python, seed 123: 1 cases, 0 differences, 0 sanitizer reports')\n"}
        r, d = self._run("engine-fuzz", "Engine generator (1,500 cases)", passing)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(d, "summary.md"), encoding="utf-8") as fh:
            self.assertIn("engine fuzz, coal-cpp vs python, seed 123", fh.read())
        r, _ = self._run("engine-fuzz", "Engine generator (1,500 cases)", {tool: "raise SystemExit(1)\n"})
        self.assertNotEqual(r.returncode, 0)
        san = "Engine generator under the sanitizers (300 cases)"
        r, _ = self._run("sanitizers", san, passing)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        r, _ = self._run("sanitizers", san, {tool: "print('x.cpp:1:1: runtime error: shift exponent')\n"})
        self.assertEqual(r.returncode, 1)
        self.assertIn("a sanitizer reported", r.stdout + r.stderr)


class TestCiJobs(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
            self.ci = fh.read()

    def _job(self, name):
        m = re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.ci, re.S | re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def test_extras_installed_and_required_in_separate_jobs(self):
        gui, mcp = self._job("gui"), self._job("mcp")
        self.assertIn('pip install -e . -e "apps/coal-app[gui]"', gui)
        self.assertIn("COAL_REQUIRE_EXTRAS: gui", gui)
        self.assertNotIn("[mcp]", gui.replace("gradio[mcp]", ""))
        self.assertIn('pip install -e . -e "apps/coal-app[mcp]"', mcp)
        self.assertIn("COAL_REQUIRE_EXTRAS: mcp", mcp)
        self.assertNotIn("[gui]", mcp)

    def test_cpp_job_cannot_go_green_without_exercising_cpp(self):
        cpp = self._job("cpp")
        # coal-cpp builds from the root CMakeLists.txt and tests/test_differential.py reads COAL_IMPL
        self.assertIn("cmake -S . -B build", cpp)
        self.assertIn("COAL_IMPL: cpp", cpp)
        self.assertIn("tests.test_differential", cpp)
        self.assertIn("tests/test_cpp_vs_python.py", cpp)
        self.assertIn("skipped", cpp)
        self.assertNotIn("not applicable", cpp.lower())

    def test_require_switch_turns_skips_into_failures(self):
        import subprocess
        import sys
        env = dict(os.environ, COAL_REQUIRE_EXTRAS="gui")
        r = subprocess.run([sys.executable, "-m", "unittest", "test_gui"], capture_output=True, text=True,
                           encoding="utf-8", cwd=os.path.join(ROOT, "apps", "coal-app", "tests"), env=env, timeout=300)
        try:
            import gradio  # noqa: F401
            has = True
        except ImportError:
            has = False
        if has:
            self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        else:
            self.assertNotEqual(r.returncode, 0)
            self.assertNotIn("skipped", r.stderr)


if __name__ == "__main__":
    unittest.main()
