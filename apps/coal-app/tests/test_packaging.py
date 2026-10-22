import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))  # the repository root


class TestPackaging(unittest.TestCase):
    def test_ci_workflow_exists(self):
        self.assertTrue(os.path.isfile(os.path.join(ROOT, ".github", "workflows", "ci.yml")))

    def test_app_imports_without_gradio(self):
        import coal_app  # noqa: F401  (must not need gradio)

    def test_sdists_leave_the_tests_out_and_the_app_carries_the_licence(self):
        # the tests need the repository (helpers, fixtures, examples, coal-cpp), so shipping them alone is useless
        for manifest in (os.path.join(ROOT, "MANIFEST.in"), os.path.join(ROOT, "apps", "coal-app", "MANIFEST.in")):
            with open(manifest, encoding="utf-8") as fh:
                self.assertIn("prune tests", fh.read(), manifest)
        with open(os.path.join(ROOT, "LICENSE"), encoding="utf-8") as fh, \
                open(os.path.join(ROOT, "apps", "coal-app", "LICENSE"), encoding="utf-8") as fh2:
            self.assertEqual(fh.read(), fh2.read(), "apps/coal-app/LICENSE must be a copy of the repository's")


def _can_build_offline():
    return subprocess.run([sys.executable, "-m", "pip", "--version"], capture_output=True,
                          timeout=60).returncode == 0


@unittest.skipUnless(_can_build_offline(), "needs pip")
class TestInstalledWheel(unittest.TestCase):
    """APP-C12: examples must be found when COAL is installed from a wheel, not only in a repo checkout."""

    def test_examples_found_after_wheel_install(self):
        t = tempfile.mkdtemp(prefix="coal_test_wheel_")
        self.addCleanup(shutil.rmtree, t, True)
        src = os.path.join(t, "src")
        shutil.copytree(ROOT, src, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.egg-info", "build",
                                                                 "PAL2NAL", "tests", "benchmarks"))
        env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "COAL_EXAMPLES_DIR")}
        # offline: COAL_WHEELHOUSE = a folder with a setuptools wheel; otherwise pip fetches the build backend
        house = os.environ.get("COAL_WHEELHOUSE")
        extra = ["--no-index", "--find-links", house] if house else []
        # two wheels since the split: the engine (coal-align, which carries the examples) and the app (coal-align-app)
        r = subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", *extra, "-w", os.path.join(t, "dist"),
                            src, os.path.join(src, "apps", "coal-app")], capture_output=True, text=True,
                           encoding="utf-8", timeout=300, env=env)
        if r.returncode != 0 and not house and re.search(r"No matching distribution|Could not find a version|"
                                                         r"Failed to establish|ProxyError|403", r.stdout + r.stderr):
            self.skipTest("no network to fetch the build backend; set COAL_WHEELHOUSE to a folder with a "
                          "setuptools wheel")
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        whls = sorted(glob.glob(os.path.join(t, "dist", "*.whl")))
        self.assertEqual(len(whls), 2, whls)
        import zipfile
        for w in whls:  # both wheels carry the MIT licence (the app wheel had none)
            names = zipfile.ZipFile(w).namelist()
            self.assertTrue(any(re.match(r"coal_align(_app)?-[^/]+\.dist-info/licenses/LICENSE$", n) for n in names),
                            (w, [n for n in names if "dist-info" in n]))
        venv = os.path.join(t, "venv")
        subprocess.run([sys.executable, "-m", "venv", "--without-pip", venv], check=True, timeout=120)
        py = os.path.join(venv, "Scripts", "python.exe") if os.name == "nt" else os.path.join(venv, "bin", "python")  # R2-D1
        r = subprocess.run([sys.executable, "-m", "pip", "--python", py, "install", "--no-deps", "--no-index", *whls],
                           capture_output=True, text=True, encoding="utf-8", timeout=300, env=env)
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        r = subprocess.run([py, "-m", "coal_app", "--check"], capture_output=True, text=True, encoding="utf-8",
                           timeout=120, env=env, cwd=t)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("examples: clean, mito, messy_ids, one_broken_gene", r.stdout)
        r = subprocess.run([py, "-c", "from coal_app import core; print(core.EXAMPLES_DIR)"], capture_output=True,
                           text=True, encoding="utf-8", timeout=60, env=env, cwd=t)
        self.assertTrue(os.path.realpath(r.stdout.strip()).startswith(os.path.realpath(venv)), r.stdout)

    def test_env_var_overrides(self):
        r = subprocess.run([sys.executable, "-c", "from coal_app import core; print(core.EXAMPLES_DIR)"],
                           capture_output=True, text=True, encoding="utf-8", timeout=60,
                           env=dict(os.environ, COAL_EXAMPLES_DIR=os.path.join(ROOT, "examples", "..", "examples")))
        self.assertEqual(os.path.realpath(r.stdout.strip()), os.path.realpath(os.path.join(ROOT, "examples")))


def _has_build():
    try:
        import build  # noqa: F401
        return True
    except ImportError:
        return False


@unittest.skipUnless(_has_build(), "needs the 'build' package (pip install build)")
class TestSdists(unittest.TestCase):
    """What `python -m build --sdist` ships for the engine and the app: no tests (they need the repository), the
    licence in both."""

    def test_sdist_contents(self):
        import tarfile
        t = tempfile.mkdtemp(prefix="coal_test_sdist_")
        self.addCleanup(shutil.rmtree, t, True)
        src = os.path.join(t, "src")
        shutil.copytree(ROOT, src, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.egg-info", "build",
                                                                 "PAL2NAL", "benchmarks", ".venv*"))
        env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH",)}
        for sub, name in (("", "coal_align"), (os.path.join("apps", "coal-app"), "coal_align_app")):
            r = subprocess.run([sys.executable, "-m", "build", "--sdist", "-o", os.path.join(t, "dist"),
                                os.path.join(src, sub)], capture_output=True, text=True, encoding="utf-8",
                               timeout=600, env=env)
            if r.returncode != 0 and re.search(r"No matching distribution|Could not find a version|"
                                               r"Failed to establish|ProxyError|403", r.stdout + r.stderr):
                self.skipTest("no network to fetch the build backend")
            self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
            [sdist] = glob.glob(os.path.join(t, "dist", name + "-*.tar.gz"))
            names = [n.split("/", 1)[1] for n in tarfile.open(sdist).getnames() if "/" in n]
            self.assertEqual([n for n in names if n.startswith("tests/")], [], (sdist, names))
            self.assertIn("LICENSE", names, (sdist, names))
            self.assertIn("pyproject.toml", names)


if __name__ == "__main__":
    unittest.main()
