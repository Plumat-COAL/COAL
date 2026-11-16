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

CPP_STEP = "C++ engine tests, every module in tests/ (must not skip for want of coal-cpp)"
SAN_STEP = "C++ tests under the sanitizers"
# The core C++ modules: no skip there but the platform ones (ci.yml), and the sanitizer job runs them the same way
CORE_CPP = ("test_differential", "test_cpp_vs_python", "test_cpp_cli", "test_cpp_readers", "test_cpp_tables",
            "test_cpp_pool")
# What the sanitizer job leaves out of ci.yml's C++ test set, and why
NOT_UNDER_SANITIZERS = {"test_wheels": "packaging (the platform wheel installs): no sanitizer business"}
# The modules a review found never run in any workflow (M1), and the ones run only by the nightly sanitizer job:
# named here so a change to how the module set is worked out cannot drop them quietly
ONCE_UNRUN = ("test_deep_findings", "test_engine_cli", "test_new_user_ux", "test_pipeline_msgs",
              "test_integration1_findings", "test_stress_findings", "test_cpp_vs_reference", "test_engine_fuzz",
              "test_orthologr_standin", "test_pwalign_matrices", "test_fuzz_harness", "test_python_api",
              "test_report_schema", "test_batch", "test_cpp_error_precedence", "test_robustness", "test_engine_ux")
# A test module that needs coal-cpp (or coal-cpp-align) for some of its tests: a "not built"-style skip, or the
# shared helpers that skip (or fail, under COAL_REQUIRE_CPP) without the binary
NEEDS_CPP = re.compile(r"not built|\bneeds_cpp\b|\bno_binary\b|\bcpp_hooks\b|COAL_CPP_BIN")


def discovered_modules() -> list[str]:
    """The modules `python -m unittest discover -s tests` collects (its default pattern, test*.py)."""
    return sorted(f[:-3] for f in os.listdir(os.path.join(ROOT, "tests")) if re.fullmatch(r"test\w*\.py", f))


def modules_needing_cpp() -> set[str]:
    out = set()
    for mod in discovered_modules():
        with open(os.path.join(ROOT, "tests", mod + ".py"), encoding="utf-8") as fh:
            if NEEDS_CPP.search(fh.read()):
                out.add(mod)
    return out


def modules_run_by(script: str) -> set[str]:
    """The test modules a step's command runs: for `discover -s tests`, every module discover collects with its
    pattern (-p / --pattern, default test*.py), less ci_cpp_tests.py's --exclude; else the tests.test_* modules it
    names. A -k filter (or --k / -k=) selects tests by name, so which modules still run cannot be told from here:
    none is counted (the guards then fail, as they should)."""
    import fnmatch
    cmd = script.replace("\\\n", " ")
    line = next((x for x in cmd.splitlines() if "ci_cpp_tests.py" in x or "-m unittest" in x), "")
    if re.search(r"(^|\s)(-k|--k)(\s|=)|\s-k\S", line):
        return set()
    excluded = {m.strip().removeprefix("tests.") for e in re.findall(r"--exclude[ =]['\"]?([\w.,]+)", line)
                for m in e.split(",") if m.strip()}
    if re.search(r"\bdiscover -s tests(\s|$)", line):
        pat = re.search(r"(?:\s-p\s*|\s--pattern[ =])['\"]?([^\s'\"]+)", line)
        pattern = pat.group(1) if pat else "test*.py"
        files = sorted(f for f in os.listdir(os.path.join(ROOT, "tests")) if fnmatch.fnmatch(f, pattern))
        return {f[:-3] for f in files if f.endswith(".py")} & set(discovered_modules()) - excluded
    return {m.split(".", 1)[1].split(".")[0] for m in re.findall(r"\btests\.test_[\w.]+", line)} - excluded


def uncommented(text: str) -> str:
    """A workflow without its comment lines (a module named only in a comment is run nowhere)."""
    return "\n".join(x for x in text.splitlines() if not x.lstrip().startswith("#"))


# Files of tests/ with a "not built" skip that discover does not collect and no workflow runs, and why
RUN_BY_HAND = {
    "fuzz_kaks_vs_reference.py": "needs KaKs_Calculator 2.0 and PAML yn00 builds that print 12+ significant digits "
                                 "(docs/kaks.md \"Validation\"), which CI does not build",
}
# coal-cpp tests that skip in ci.yml's C++ job for want of something other than the binary, and the nightly job
# that runs them with it (nightly.yml cpp-extras): module -> what it needs
RUN_NIGHTLY = {
    "test_cpp_vs_reference": "PAL2NAL (TestCppVsPal2nalFuzz) and COAL_CPP_NOHOOKS (the hookless-build test)",
    "test_stress_findings": "COAL_SLOW (QuadraticExactFitScan, UnrelatedTwentyThousand)",
    "test_align": "Biopython (LiveBiopython)",
    "test_robustness": "root, for a 1 MB tmpfs (TestDiskFull)",
}
# coal-cpp tests that run in no job at all, and why (keep this short: each is a gap)
NEVER_RUN = {
    "test_kaks_myn_kappa_inf: the C++ API part on Windows": "the MSVC build has coal_core.lib, not libcoal_core.a "
                                                          "(it runs on Linux and macOS)",
}


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
        self.assertEqual(sorted(jobs), ["cpp-extras", "engine-fuzz", "kaks-properties", "msa-tools", "pal2nal-fuzz",
                                        "sanitizers", "tsan", "tutorials"])
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

    def test_msa_tools_job(self):
        """ci.yml installs Clustal Omega and MAFFT best effort, so its MSA tests may skip: here the install must
        succeed and those tests must not skip for want of the tools (REVIEW-COAL-3 L4)."""
        job = self._job("msa-tools")
        self.assertIn("python -m coal pipeline install clustalo mafft\n", job)
        self.assertIn("assert c and m", job)
        self.assertNotIn("|| {", job.split("Install Clustal Omega", 1)[1].split("- name:", 1)[0])   # no fallback
        self.assertIn('python .github/ci_cpp_tests.py --no-skip "$mods"', job)
        mods = re.search(r"mods=(\S+)", job).group(1).split(",")
        for m in ("test_msa", "test_dnds_pipeline", "test_pipeline_provenance"):
            self.assertIn(m, mods)
        for m in mods:
            self.assertTrue(os.path.isfile(os.path.join(ROOT, "tests", m + ".py")), m)
        allow = re.search(r"--allow '([^']+)'", job).group(1)
        # every skip reason may stay, except one that names the MSA tools
        rx = re.compile(allow, re.I)
        for reason in ("needs coal-cpp and Clustal Omega", "needs MAFFT", "neither clustalo nor mafft found"):
            self.assertIsNone(rx.search(reason), reason)
        for reason in ("POSIX shell stand-ins", "symbolic links need a privilege on Windows", "Biopython missing"):
            self.assertIsNotNone(rx.search(reason), reason)
        for env in ("COAL_REQUIRE_CPP", "COAL_REQUIRE_ALIGN", "COAL_REQUIRE_DIAMOND=1"):
            self.assertIn(env, job)
        # ci.yml's install stays best effort (it says so) and this job is what makes up for it
        ci_step = re.search(r"\n((?:      #[^\n]*\n)*)      - name: Install Clustal Omega and MAFFT.*?\n(?=      - name:)",
                            self.ci, re.S).group(0)
        self.assertIn("continue-on-error: true", ci_step)
        self.assertIn("nightly.yml msa-tools", ci_step)

    def test_sanitizer_job(self):
        job = self._job("sanitizers")
        self.assertIn("-fsanitize=address,undefined -fno-omit-frame-pointer -O1 -g", job)
        self.assertIn('-DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"', job)
        self.assertIn("-DCOAL_WARNINGS_AS_ERRORS=ON", job)
        self.assertIn("ASAN_OPTIONS: detect_leaks=1", job)
        self.assertIn("UBSAN_OPTIONS: print_stacktrace=1:halt_on_error=1", job)
        self.assertIn('COAL_SKIP_MEMORY_TESTS: "1"', job)
        self.assertIn("COAL_IMPL: cpp", job)
        san = modules_run_by(step_script("sanitizers", SAN_STEP, "nightly.yml"))
        for mod in CORE_CPP + ("test_stdin", "test_robustness", "test_cpp_error_precedence", "test_batch"):
            self.assertIn(mod, san, mod)
        self.assertIn("--no-skip " + ",".join(CORE_CPP) + " ", job)
        # DIAMOND installed as in ci.yml's C++ job, before the tests (its tests must not skip there either)
        install = re.search(r"(?m)^        run: python -c .*install_diamond.*$", job).group(0)
        self.assertIn(install, _ci())
        self.assertLess(job.index(install), job.index(f"- name: {SAN_STEP}"))
        self.assertIn('engine_fuzz.py --n 300 --seed "$FUZZ_SEED" --sanitizer', job)
        self.assertIn('tests/test_kaks_properties.py --cpp build-asan/coal-cpp --n 300 --seed "$FUZZ_SEED"', job)
        self.assertEqual(job.count('grep -qE "AddressSanitizer|LeakSanitizer|runtime error"'), 3)

    def test_kaks_properties_job(self):
        """The Ka/Ks property fuzz with tonight's seed: COAL_SLOW unit tests, then the script's campaign, which
        fails on any finding that is not a known bug; the findings are uploaded on failure."""
        job = self._job("kaks-properties")
        self.assertIn('COAL_SLOW: "1"', job)
        self.assertIn('COAL_REQUIRE_CPP: "1"', job)
        self.assertIn('KAKS_FUZZ_SEED="$FUZZ_SEED" python -m unittest -v tests.test_kaks_properties', job)
        self.assertIn('tests/test_kaks_properties.py --n 1500 --seed "$FUZZ_SEED" --keep kaks-props-keep', job)
        self.assertIn("if: failure()", job)
        self.assertIn("kaks-props-keep/", job)
        self.assertTrue(os.path.isfile(os.path.join(ROOT, "tests", "test_kaks_properties.py")))

    def test_tsan_job(self):
        job = self._job("tsan")
        self.assertIn('-DCMAKE_CXX_FLAGS="-fsanitize=thread -g -O1 -Wno-error=free-nonheap-object"', job)
        # GCC 15's false positive under -fsanitize=thread: waived in this job only, and only that one warning
        everywhere = "\n".join(x for x in (self.wf + self.ci).splitlines() if not x.lstrip().startswith("#"))
        self.assertEqual(everywhere.count("-Wno-error="), 1)
        self.assertNotRegex(everywhere, r"-Wno-(?!error=free-nonheap-object)[\w-]+")
        self.assertIn('-DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread"', job)
        self.assertIn("-DCOAL_WARNINGS_AS_ERRORS=ON", job)
        self.assertIn("test -x build-tsan/coal-cpp-align", job)
        self.assertIn("nm build-tsan/coal-cpp | grep -q __tsan_init", job)       # a positive control
        self.assertIn("nm build-tsan/coal-cpp-align | grep -q __tsan_init", job)
        self.assertIn('tests/tools/tsan_threads.py --build build-tsan --seed "$FUZZ_SEED" --jobs 4', job)
        self.assertIn('grep -q "ThreadSanitizer" tsan.txt', job)
        self.assertIn("if: failure()", job)
        self.assertIn("tsan-keep/", job)
        # the tool runs kaks, dnds and codon-dnds with all nine methods and both flavours, and the aligner
        with open(os.path.join(ROOT, "tests", "tools", "tsan_threads.py"), encoding="utf-8") as fh:
            tool = fh.read()
        for word in ('"Comeron", "Li"', '"kaks"', '"dnds"', '"codon-dnds"', "coal-cpp-align", '"gestimator"',
                     '"seqinr"', '"@J"'):
            self.assertIn(word, tool, word)

    def test_tutorials_job_is_cis_without_the_paths_filter(self):
        """The nightly tutorials job runs ci.yml's tutorials steps (the same scripts, step for step) every night,
        without the first step that skips the job when nothing under the tutorials' paths changed."""
        job = self._job("tutorials")
        self.assertNotIn("steps.changes", job)
        self.assertNotIn("if: github.event_name", job)
        steps = ("Install the engine, build coal-cpp and coal-cpp-align", "Install DIAMOND, Clustal Omega, MAFFT and PAL2NAL v14",
                 "Download the tutorial data (NCBI, SHA-256 checked)", "Run the tutorials and compare with what they show")
        for step in steps:
            self.assertEqual(step_script("tutorials", step, "nightly.yml"), step_script("tutorials", step), step)
        self.assertIn("scripts/run_tutorials.py --full", job)
        self.assertIn("key: tutorial-data-${{ hashFiles('docs/tutorials/data/*.sh') }}", job)
        self.assertIn("if: failure()", job)
        self.assertIn("retention-days:", job)

    def test_the_ci_cpp_test_set_is_the_nightly_one(self):
        """The sanitizer job runs ci.yml's C++ test set (every module in tests/), less NOT_UNDER_SANITIZERS: every
        module that needs coal-cpp is tested under ASan and UBSan too."""
        ci_mods = modules_run_by(step_script("cpp", CPP_STEP))
        san_mods = modules_run_by(step_script("sanitizers", SAN_STEP, "nightly.yml"))
        self.assertEqual(ci_mods - san_mods, set(NOT_UNDER_SANITIZERS))
        self.assertLessEqual(modules_needing_cpp() - set(NOT_UNDER_SANITIZERS), san_mods)
        for mod in NOT_UNDER_SANITIZERS:
            self.assertIn(mod, ci_mods)


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
        with open(os.path.join(ROOT, ".github", "ci_cpp_tests.py"), encoding="utf-8") as fh:
            files = {".github/ci_cpp_tests.py": fh.read()}
        for mod in CORE_CPP + ("test_batch", "test_wheels"):
            files[f"tests/{mod}.py"] = "import unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n        pass\n"
        r, _ = self._run("sanitizers", SAN_STEP, files)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("(not run: test_wheels)", r.stdout + r.stderr)
        leaky = dict(files)
        leaky["tests/test_batch.py"] = ("import sys, unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n"
                                        "        print('==7==ERROR: LeakSanitizer: detected memory leaks', file=sys.stderr)\n")
        r, _ = self._run("sanitizers", SAN_STEP, leaky)
        self.assertEqual(r.returncode, 1)
        self.assertIn("a sanitizer reported", r.stdout + r.stderr)
        unbuilt = dict(files)
        unbuilt["tests/test_batch.py"] = ("import unittest\nclass T(unittest.TestCase):\n"
                                          "    @unittest.skip('coal-cpp not built')\n    def test_x(self):\n        pass\n")
        r, _ = self._run("sanitizers", SAN_STEP, unbuilt)
        self.assertEqual(r.returncode, 1)
        self.assertIn("not built", r.stdout + r.stderr)

    def test_kaks_property_steps_fail_on_a_finding_and_write_the_summary(self):
        tool = "tests/test_kaks_properties.py"

        def stand_in(script_status, unit_ok=True, extra=""):
            return ("import sys, unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n"
                    f"        self.assertTrue({unit_ok})\n"
                    "if __name__ == '__main__':\n"
                    f"    {extra}print('kaks properties: 10 rows, 10 swap; findings none (1 s, seed 123, 1 pairs)')\n"
                    f"    sys.exit({script_status})\n")
        step = "Ka/Ks properties (COAL_SLOW unit tests and a 1,500-pair campaign)"
        files = {"tests/__init__.py": ""}
        r, d = self._run("kaks-properties", step, {**files, tool: stand_in(0)})
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(d, "summary.md"), encoding="utf-8") as fh:
            self.assertIn("kaks properties: 10 rows", fh.read())
        r, _ = self._run("kaks-properties", step, {**files, tool: stand_in(1)})          # a finding of the campaign
        self.assertNotEqual(r.returncode, 0)
        r, _ = self._run("kaks-properties", step, {**files, tool: stand_in(0, unit_ok=False)})   # a failed unit test
        self.assertNotEqual(r.returncode, 0)
        san = "Ka/Ks property harness under the sanitizers (300 pairs)"
        r, _ = self._run("sanitizers", san, {**files, tool: stand_in(0)})
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        leak = "print('==7==ERROR: LeakSanitizer: detected memory leaks', file=sys.stderr); "
        r, _ = self._run("sanitizers", san, {**files, tool: stand_in(0, extra=leak)})
        self.assertEqual(r.returncode, 1)
        self.assertIn("a sanitizer reported", r.stdout + r.stderr)
        r, _ = self._run("sanitizers", san, {**files, tool: stand_in(1)})
        self.assertNotEqual(r.returncode, 0)

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


@unittest.skipUnless(BASH and os.name != "nt", "needs a POSIX bash (the step runs on ubuntu-latest)")
class TestCppTestStep(unittest.TestCase):
    """ci.yml's C++ test step, run for real with stand-in test modules: a "not built" skip anywhere fails, any skip
    in a core C++ module fails unless it is a platform one, a core module that is gone fails, and a skip that
    a test's own output surrounds is still seen."""

    OK = "import unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n        pass\n"
    # the Windows job's matrix.allow (the widest of the three)
    allow = re.search(r"- os: windows-latest\n\s+allow: \"([^\"]*)\"", _ci()).group(1)

    @staticmethod
    def skipping(reason, noisy=False):
        say = "        import sys; print('chatter', file=sys.stderr)\n" if noisy else ""
        return ("import unittest\nclass T(unittest.TestCase):\n    def test_x(self):\n" + say +
                f"        self.skipTest({reason!r})\n")

    def run_step(self, files):
        d = self.dir = tempfile.mkdtemp(prefix="coal_ci_step_")
        self.addCleanup(shutil.rmtree, d, True)
        with open(os.path.join(ROOT, ".github", "ci_cpp_tests.py"), encoding="utf-8") as fh:
            body = {".github/ci_cpp_tests.py": fh.read()}
        body.update({f"tests/{m}.py": self.OK for m in CORE_CPP + ("test_wheels",)})
        body.update(files)
        for name, text in body.items():
            os.makedirs(os.path.dirname(os.path.join(d, name)), exist_ok=True)
            with open(os.path.join(d, name), "w", encoding="utf-8") as fh:
                fh.write(text)
        env = dict(os.environ, GITHUB_STEP_SUMMARY=os.path.join(d, "summary.md"), CPP_ALLOW=self.allow,
                   CPP_NO_SKIP_EXTRA="")
        script = re.sub(r"(?m)^python ", __import__("sys").executable + " ", step_script("cpp", CPP_STEP))
        r = subprocess.run([BASH, "--noprofile", "--norc", "-eo", "pipefail", "-c", script], cwd=d, env=env,
                           capture_output=True, text=True, timeout=120)
        return r, r.stdout + r.stderr

    def test_passing_run_passes_and_writes_the_summary(self):
        r, out = self.run_step({"tests/test_other.py": self.OK})
        self.assertEqual(r.returncode, 0, out)
        with open(os.path.join(self.dir, "summary.md"), encoding="utf-8") as fh:
            self.assertIn("C++ tests: 8 run in 8 modules, 0 skipped", fh.read())

    def test_not_built_skip_in_any_module_fails(self):
        for reason in ("coal-cpp not built", "build/coal-cpp (COAL_CPP_BIN) is not built",
                       "coal-cpp-align not built (set COAL_CPP_ALIGN)", "diamond or coal-cpp-align not available"):
            with self.subTest(reason=reason):
                r, out = self.run_step({"tests/test_pipeline_x.py": self.skipping(reason, noisy=True)})
                self.assertEqual(r.returncode, 1, out)
                self.assertIn("skipped for want of a binary: test_pipeline_x.T.test_x", out)

    def test_other_skips_outside_the_core_modules_pass(self):
        r, out = self.run_step({"tests/test_pipeline_x.py": self.skipping("needs coal-cpp and Clustal Omega")})
        self.assertEqual(r.returncode, 0, out)

    def test_core_module_skips_fail_but_the_platform_ones(self):
        for reason, rc in (("POSIX signals", 0), ("POSIX hard links and signals", 0),
                           ("peak memory read with getrusage (Linux: KB)", 0),
                           ("live comparison needs Python 3.13's argparse", 1), ("POSIX shell stand-in", 1),
                           ("set PAL2NAL=/path/to/pal2nal.pl", 1)):
            for mod in ("test_cpp_cli", "test_wheels"):
                with self.subTest(reason=reason, mod=mod):
                    r, out = self.run_step({f"tests/{mod}.py": self.skipping(reason)})
                    self.assertEqual(r.returncode, rc, out)
                    if rc:
                        self.assertIn(f"skipped in a module that must not skip ({mod})", out)

    def test_a_core_module_that_is_gone_fails(self):
        r, out = self.run_step({"tests/test_cpp_pool.py": "# no tests\n"})
        self.assertEqual(r.returncode, 1, out)
        self.assertIn("test_cpp_pool: no test of this module ran", out)

    def test_a_run_of_no_tests_fails(self):
        files = {f"tests/{m}.py": "# no tests\n" for m in CORE_CPP + ("test_wheels",)}
        r, out = self.run_step(files)
        self.assertEqual(r.returncode, 1, out)
        self.assertIn("no test ran at all", out)

    def test_a_failing_test_fails_with_an_annotation(self):
        r, out = self.run_step({"tests/test_pipeline_x.py": self.OK.replace("pass", "self.fail('boom')")})
        self.assertEqual(r.returncode, 1, out)
        self.assertIn("::error title=FAIL: test_pipeline_x.T.test_x::", out)


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

    def test_every_module_that_needs_coal_cpp_runs_in_the_cpp_job(self):
        """M1: the C++ job (all three systems) runs every module of tests/ with a coal-cpp skip, so no such test
        is skipped everywhere (the Python jobs build no binary). It runs them with the binary required
        (COAL_REQUIRE_CPP, COAL_REQUIRE_ALIGN; COAL_REQUIRE_DIAMOND from the DIAMOND step, which comes first and
        has no `if:`), through .github/ci_cpp_tests.py (a "not built" skip fails)."""
        cpp = self._job("cpp")
        self.assertIn("os: [ubuntu-latest, windows-latest, macos-latest]", cpp)
        step = re.search(rf"      - name: {re.escape(CPP_STEP)}\n(.*?)(?=\n      - name:|\n      # The release|\Z)",
                         cpp, re.S)
        self.assertIsNotNone(step)
        self.assertNotRegex(step.group(1), r"(?m)^        if:")
        for env in ('COAL_REQUIRE_CPP: "1"', 'COAL_REQUIRE_ALIGN: "1"', "COAL_IMPL: cpp"):
            self.assertIn(env, step.group(1))
        script = step_script("cpp", CPP_STEP)
        self.assertRegex(script, r"(?m)^python \.github/ci_cpp_tests\.py ")
        run = modules_run_by(script)
        needing = modules_needing_cpp()
        self.assertGreater(len(needing), 40)          # the detection still finds them
        self.assertEqual(needing - run, set(), "modules with coal-cpp skips that the C++ job does not run")
        for mod in ONCE_UNRUN + ("test_kaks", "test_kaks_comeron_li", "test_align", "test_msa", "test_dnds_pipeline",
                                 "test_short_options", "test_output_lock", "test_kaks_dispatch", "test_wheels"):
            self.assertIn(mod, run, mod)
            self.assertIn(mod, needing | {"test_msa"}, mod)
            self.assertTrue(os.path.isfile(os.path.join(ROOT, "tests", mod + ".py")), mod)
        self.assertIn('--no-skip "' + ",".join(CORE_CPP) + ',test_wheels$CPP_NO_SKIP_EXTRA" ', script.replace("\\\n", " "))
        self.assertIn('--allow "$CPP_ALLOW"', script)
        order = [m.group(1) for m in re.finditer(r"(?m)^      - name: (.*)$", cpp)]
        diamond = next(i for i, n in enumerate(order) if n.startswith("Install DIAMOND (pinned release"))
        self.assertLess(diamond, order.index(CPP_STEP))     # COAL_REQUIRE_DIAMOND reaches the tests
        self.assertLess(next(i for i, n in enumerate(order) if n.startswith("Install Clustal Omega and MAFFT")),
                        order.index(CPP_STEP))
        d = re.search(r"      - name: Install DIAMOND \(pinned release.*?\n(.*?)(?=\n      (- name:|#))", cpp, re.S)
        self.assertIsNotNone(d)
        self.assertNotIn("if:", d.group(1))
        self.assertIn("COAL_REQUIRE_DIAMOND=1", d.group(1))
        self.assertIn("--align-binary build-manylinux/coal-cpp-align", cpp)
        self.assertIn("packaging/wheels/check_wheel_binaries.py", cpp)

    def test_allowed_skips_are_per_system(self):
        """The skips allowed in the core modules are each system's own (matrix.allow): none on Linux, the
        Linux-only peak-memory test on macOS, that and POSIX signals / pipes / hard links on Windows; nothing for
        an argparse that differs (the job runs Python 3.13, where the live comparison runs). On Linux the dispatch
        test is a core module too, and the AVX2 copy must be the one that runs."""
        cpp = self._job("cpp")
        rows = dict(re.findall(r"- os: ([\w-]+)\n\s+allow: \"([^\"]*)\"", cpp))
        self.assertEqual(rows, {"ubuntu-latest": "", "macos-latest": "peak memory read with getrusage",
                                "windows-latest": "peak memory read with getrusage|POSIX (signals|pipes|hard links)"})
        self.assertRegex(cpp, r"- os: ubuntu-latest\n\s+allow: \"\"\n\s+no_skip_extra: \",test_kaks_dispatch\"")
        self.assertEqual(len(re.findall(r"no_skip_extra: \"\"", cpp)), 2)
        self.assertNotIn("live comparison", uncommented(_ci() + _ci("nightly.yml")))
        self.assertIn('test "$("$b" dnds --kaks-isa)" = avx2', step_script("cpp", "Point the tests at the binary"))
        self.assertIn('--allow "peak memory read with getrusage"', step_script("sanitizers", SAN_STEP, "nightly.yml"))

    def test_the_nightly_runs_the_cpp_tests_ci_cannot(self):
        """RUN_NIGHTLY: the coal-cpp tests that skip in the C++ job for want of PAL2NAL, a hookless build, COAL_SLOW,
        Biopython or root run in nightly.yml's cpp-extras job, with what they need, and must not skip there."""
        wf = _ci("nightly.yml")
        job = re.search(r"^  cpp-extras:\n(.*?)(?=^  [\w-]+:\n|\Z)", wf, re.S | re.M).group(1)
        main = step_script("cpp-extras", "PAL2NAL fuzz, hookless build, COAL_SLOW and Biopython tests (must not skip)",
                           "nightly.yml")
        root = step_script("cpp-extras", "Full file system tests as root (must not skip)", "nightly.yml")
        run = modules_run_by(main) | modules_run_by(root)
        self.assertEqual(set(RUN_NIGHTLY) - run, set())
        for mod in RUN_NIGHTLY:
            self.assertRegex(main + root, rf"--no-skip [\w,]*{mod}\b", mod)
        for need in ('COAL_SLOW: "1"', "PAL2NAL: PAL2NAL/pal2nal.v14/pal2nal.pl", "COAL_CPP_NOHOOKS=",
                     "pip install -e . biopython psutil", "run: python scripts/fetch_pal2nal.py",
                     'COAL_REQUIRE_CPP: "1"', "COAL_REQUIRE_DIAMOND=1", "-DCOAL_TEST_HOOKS=ON"):
            self.assertIn(need, job, need)
        # sudo-rs ignores --preserve-env: what the tests need is passed by name, and TestDiskFull fails without coal-cpp
        self.assertNotRegex(root, r"(?m)^\s*sudo\b.*--preserve-env")
        self.assertIn('sudo env "PATH=$PATH" "COAL_CPP_BIN=$COAL_CPP_BIN" COAL_REQUIRE_CPP=1 COAL_IMPL=cpp', root)
        self.assertIn("tests.test_robustness.TestDiskFull", root)
        self.assertTrue(NEVER_RUN)   # the gaps that stay are named, with their reasons
        for what, why in NEVER_RUN.items():
            self.assertTrue(why, what)

    def test_sanitizer_job_has_psutil_and_the_msa_tools(self):
        job = re.search(r"^  sanitizers:\n(.*?)(?=^  [\w-]+:\n|\Z)", _ci("nightly.yml"), re.S | re.M).group(1)
        self.assertIn("run: python -m pip install -e . psutil", job)
        steps = [m.group(1) for m in re.finditer(r"(?m)^      - name: (.*)$", job)]
        self.assertLess(steps.index("Install Clustal Omega and MAFFT (Bioconda, by coal)"), steps.index(SAN_STEP))
        self.assertIn("run: python -m coal pipeline install clustalo mafft", job)

    def test_engine_ux_runs_in_the_cpp_job(self):
        """tests.test_engine_ux (the C++ engine's messages) needs coal-cpp: it runs in the C++ job and in the
        sanitizer job."""
        self.assertIn("test_engine_ux", modules_needing_cpp())
        self.assertIn("test_engine_ux", modules_run_by(step_script("cpp", CPP_STEP)))
        self.assertIn("test_engine_ux", modules_run_by(step_script("sanitizers", SAN_STEP, "nightly.yml")))

    def test_x86_macos_job_runs_the_avx2_copy(self):
        """AppleClang on x86_64 is the only place the AVX2 Ka/Ks copy is compiled that no other job runs: the
        x86_64 macOS job builds there and the dispatch test must not skip."""
        job = self._job("cpp-macos-x86")
        self.assertIn("runs-on: macos-15-intel", job)
        self.assertIn('test "$(./build/coal-cpp dnds --kaks-isa)" = avx2', job)
        script = step_script("cpp-macos-x86", "Dispatch, Ka/Ks and golden tests with coal-cpp (must not skip)")
        self.assertIn("--no-skip test_kaks_dispatch,", script)
        self.assertIn("test_kaks_dispatch", modules_run_by(script))

    def test_files_with_cpp_skips_outside_discover_run_somewhere(self):
        """A file of tests/ with a "not built" skip that discover does not collect (the fuzz scripts) is run by
        some workflow."""
        wf = uncommented(_ci() + "\n" + _ci("nightly.yml"))
        seen = set()
        for f in sorted(os.listdir(os.path.join(ROOT, "tests"))):
            if not f.endswith(".py") or re.fullmatch(r"test\w*\.py", f):
                continue
            with open(os.path.join(ROOT, "tests", f), encoding="utf-8") as fh:
                if "not built" not in fh.read():
                    continue
            seen.add(f)
            if f in RUN_BY_HAND:
                self.assertNotIn(f"tests/{f}", wf, f"{f} runs in a workflow now: take it off RUN_BY_HAND")
            else:
                self.assertIn(f"tests/{f}", wf, f"{f}: run by no workflow (or add it to RUN_BY_HAND with a reason)")
        self.assertLessEqual(set(RUN_BY_HAND), seen)

    def test_test_hooks_only_in_test_builds(self):
        """coal-cpp's COAL_TEST_* hooks are a CMake option (default OFF): on in the builds the tests run against,
        nowhere else (release.yml, the wheels' build script, the Bioconda recipe, the app's binary)."""
        def text(*parts):
            with open(os.path.join(ROOT, *parts), encoding="utf-8") as fh:
                return fh.read()
        cmake = text("CMakeLists.txt")
        self.assertRegex(cmake, r'option\(COAL_TEST_HOOKS "[^"]+" OFF\)')
        self.assertIn("target_compile_definitions(coal-cpp PRIVATE COAL_TEST_HOOKS)", cmake)
        cpp = self._job("cpp")
        self.assertIn("cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOAL_WARNINGS_AS_ERRORS=ON -DCOAL_TEST_HOOKS=ON",
                      cpp)
        self.assertEqual(_ci("ci.yml").count("COAL_TEST_HOOKS=ON"), 2)      # that line and the comment above it
        self.assertNotIn("COAL_TEST_HOOKS", self._job("windows-app"))
        self.assertIn("-B build-asan -DCMAKE_BUILD_TYPE=Debug -DCOAL_WARNINGS_AS_ERRORS=ON -DCOAL_TEST_HOOKS=ON",
                      _ci("nightly.yml"))
        for rel in (".github/workflows/release.yml", "packaging/wheels/build_linux_binary.sh",
                    "packaging/bioconda/build.sh"):
            self.assertNotIn("COAL_TEST_HOOKS", text(*rel.split("/")), rel)
        # no hook is read outside an #ifdef COAL_TEST_HOOKS block, and only through test_hook_env
        for name in ("main.cpp", "sysio.cpp"):
            depth, src = 0, text("cpp", "src", name)
            for line in src.splitlines():
                code = line.split("//")[0]
                if code.startswith("#ifdef COAL_TEST_HOOKS"):
                    depth = 1
                elif depth and code.startswith("#if"):
                    depth += 1
                elif depth and code.startswith("#endif"):
                    depth -= 1
                elif depth == 1 and code.startswith("#else"):
                    depth = 0
                elif "COAL_TEST_" in code or "test_hook_env" in code:
                    self.assertTrue(depth, f"{name}: {line.strip()}")
        docs = text("docs", "development.md")
        for var in sorted(set(re.findall(r"COAL_TEST_[A-Z_]+", text("cpp", "src", "sysio.cpp") + text("cpp", "src", "main.cpp")))):
            self.assertIn(var, docs, f"{var} is not in docs/development.md \"Test hooks\"")

    def test_cpp_job_cannot_go_green_without_exercising_cpp(self):
        cpp = self._job("cpp")
        # coal-cpp builds from the root CMakeLists.txt and tests/test_differential.py reads COAL_IMPL
        self.assertIn("cmake -S . -B build", cpp)
        self.assertIn("COAL_IMPL: cpp", cpp)
        self.assertIn("test_differential", modules_run_by(step_script("cpp", CPP_STEP)))
        self.assertIn("tests/test_cpp_vs_python.py", cpp)
        self.assertIn("python .github/ci_cpp_tests.py", cpp)
        self.assertIn('--no-skip "test_differential,', cpp)
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


class TestRPackageWorkflow(unittest.TestCase):
    """r-package.yml: R CMD check of packaging/r/coalr on three systems and the vendored-sources check, in a
    workflow of its own with a paths filter (the main matrix is not slowed), with ci.yml's hygiene."""

    def setUp(self):
        self.wf = _ci("r-package.yml")

    def _job(self, name):
        m = re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.wf, re.S | re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def test_triggers_and_hygiene(self):
        self.assertRegex(self.wf, r"(?m)^permissions:\n  contents: read")
        self.assertIn("workflow_dispatch:", self.wf)
        head = self.wf.split("\njobs:\n", 1)[0]
        for path in ("packaging/r/**", "cpp/src/**", "scripts/sync_r_sources.py", ".github/workflows/r-package.yml"):
            self.assertEqual(head.count(f'- "{path}"'), 2, path)        # for push and for pull_request
        self.assertNotRegex(head, r"[&*]\w+\s*$")                       # no YAML anchors
        jobs = re.findall(r"^  ([\w-]+):\n", self.wf.split("\njobs:\n", 1)[1], re.M)
        self.assertEqual(sorted(jobs), ["check", "vendored-sources"])
        for j in jobs:
            self.assertRegex(self._job(j), r"(?m)^    timeout-minutes: \d+", j)
            self.assertNotIn("continue-on-error", self._job(j), j)
        ci = {u.split("@")[0]: u for u in re.findall(r"uses: (\S+)", _ci())}
        uses = re.findall(r"uses: (\S+)", self.wf)
        r_lib = set()
        for u in uses:
            self.assertRegex(u, r"@[0-9a-f]{40}$", f"unpinned action {u}")
            if u.startswith("r-lib/actions/"):
                r_lib.add(u.split("@")[1])
            else:
                self.assertEqual(u, ci.get(u.split("@")[0]), f"{u}: not the SHA ci.yml uses")
        self.assertEqual(len(r_lib), 1, "the three r-lib actions come from one commit")
        # the R jobs are not in ci.yml: the main matrix does not wait for R
        self.assertNotIn("r-lib/actions", _ci())

    def test_check_job(self):
        job = self._job("check")
        for a in ("setup-r@", "setup-r-dependencies@", "check-r-package@"):
            self.assertIn(f"r-lib/actions/{a}", job)
        self.assertEqual(job.count("working-directory: packaging/r/coalr"), 2)
        rows = re.findall(r"- \{os: ([\w.-]+), r: ([\w.-]+)\}", job)
        self.assertEqual(sorted(rows), [("macos-latest", "release"), ("ubuntu-latest", "oldrel-1"),
                                        ("ubuntu-latest", "release"), ("windows-latest", "release")])
        self.assertIn("r-version: ${{ matrix.r }}", job)
        # r-lib/actions v2.14.0: v2.9.0 failed on ubuntu-latest (pak installed through a sudo that ignores -E) and
        # hung on windows-latest (in the Rtools45 installer, which v2.9.0 predates), PR #55
        self.assertEqual(job.count("@f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05  # v2.14.0"), 4)   # with setup-pandoc
        self.assertIn("r-lib/actions/setup-pandoc@", job)                          # (the vignette needs pandoc)
        self.assertIn("    needs: vendored-sources\n", job)
        setup = job.split("r-lib/actions/setup-r@", 1)[1].split("- uses:", 1)[0]
        self.assertRegex(setup, r"timeout-minutes: \d+\n")       # an installer that waits for a click: minutes
        self.assertIn('"--as-cran"', job)
        self.assertIn('"--no-manual"', job)
        self.assertIn("""error-on: '"warning"'""", job)
        self.assertTrue(os.path.isfile(os.path.join(ROOT, "packaging", "r", "coalr", "DESCRIPTION")))

    def test_vendored_sources_job_runs_the_check(self):
        self.assertIn("run: python scripts/sync_r_sources.py --check", self._job("vendored-sources"))
        self.assertTrue(os.path.isfile(os.path.join(ROOT, "scripts", "sync_r_sources.py")))


class TestReleaseSmokeWorkflow(unittest.TestCase):
    """release-smoke.yml installs the artifacts of a release.yml run and runs them: by hand only, a read-only token
    (plus actions: read for another run's artifacts), nothing that publishes, actions at release.yml's SHAs."""

    def test_hygiene(self):
        wf = _ci("release-smoke.yml")
        self.assertRegex(wf, r"(?m)^on:\n  workflow_dispatch:\n")
        self.assertNotRegex(wf, r"(?m)^  (push|pull_request|schedule|workflow_run):")
        self.assertRegex(wf, r"(?m)^permissions:\n  contents: read\n  actions: read[^\n]*\n\n")
        self.assertNotRegex(wf, r"(?m)^[ ]+permissions:")          # no job widens it
        for word in ("gh release", "pypi-publish", "contents: write", "id-token", "secrets."):
            self.assertNotIn(word, wf)
        rel = {r.split("@")[0]: r for r in re.findall(r"(?m)^\s+(?:- )?uses: (\S+)", _ci("release.yml"))}
        uses = re.findall(r"(?m)^\s+(?:- )?uses: (\S+)", wf)
        self.assertTrue(uses)
        for u in uses:
            self.assertRegex(u, r"@[0-9a-f]{40}$", f"unpinned action {u}")
            self.assertEqual(u, rel.get(u.split("@")[0], u), "not the SHA release.yml uses")
        self.assertEqual(len(re.findall(r"timeout-minutes: \d+", wf)), 4)   # every job
        self.assertNotIn("actions/checkout@", wf.split("\n  smoke:\n", 1)[1].split("\n  compare:\n", 1)[0])


if __name__ == "__main__":
    unittest.main()
