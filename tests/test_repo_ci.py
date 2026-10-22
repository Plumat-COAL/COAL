"""The CI workflow: step scripts run for real (R2-D4 / R2-D11) instead of grepped, required jobs, pinned actions, minimal permissions, and the extras jobs."""
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
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._job("explainer"), self.re.M))

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
