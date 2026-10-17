"""WEB-1..WEB-6: the download website (site/) and the release/pages workflows behind it."""
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE = os.path.join(ROOT, "site")
WORKFLOWS = os.path.join(ROOT, ".github", "workflows")
REPO_URL = "https://github.com/plumat-coal/COAL"

try:
    import yaml  # PyYAML: comes with gradio (the gui extra); CI's gui job sets COAL_REQUIRE_EXTRAS=gui
except ImportError:  # pragma: no cover - depends on the environment
    yaml = None
REQUIRE_YAML = "gui" in os.environ.get("COAL_REQUIRE_EXTRAS", "").split(",")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "scripts", name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as fh:
        return fh.read()


def _workflow(name):
    if yaml is None:
        if REQUIRE_YAML:
            raise AssertionError("PyYAML is missing although COAL_REQUIRE_EXTRAS=gui")
        raise unittest.SkipTest("PyYAML not installed (it comes with the gui extra)")
    with open(os.path.join(WORKFLOWS, name), encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _steps(wf, job):
    return {s.get("name", s.get("uses", "")): s for s in wf["jobs"][job]["steps"]}


def release_file_names(version="9.9.9"):
    """The exact file list the release job checks (step "Check the release files"), with ${V} filled in."""
    run = _steps(_workflow("release.yml"), "release")["Check the release files"]["run"]
    block = re.search(r"expected=\((.*?)\)", run, re.S).group(1)
    return [n.replace("${V}", version) for n in re.findall(r'"([^"]+)"', block)]


class TestReleaseScripts(unittest.TestCase):
    def test_version_check(self):
        rv = _load("release_version")
        init, proj = rv.source_versions()
        self.assertEqual(init, proj, "src/coal/__init__.py and pyproject.toml disagree")
        self.assertEqual(rv.check("v" + init), (init, False))
        for bad in ("v9.9.9", init, "v" + init + "x", "", "release-" + init):
            with self.assertRaises(ValueError, msg=bad):
                rv.check(bad)
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "release_version.py"), "v9.9.9"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 1)
        self.assertIn("bump both", r.stderr)

    def test_prerelease_tags(self):
        rv = _load("release_version")
        self.assertTrue(rv.TAG_RE.match("v1.0.0rc1"))
        self.assertFalse(rv.TAG_RE.match("v1.0.0-rc1"))

    def test_notes(self):
        rn = _load("release_notes")
        notes = rn.fixed_part("1.2.3")
        for must in ("not code-signed", "Windows protected your PC", "Smart App Control", "Run anyway",
                     "pip install", "coal-1.2.3-py3-none-any.whl", "releases/download/v1.2.3/", "SHA256SUMS.txt",
                     "Get-FileHash", "sha256sum"):
            self.assertIn(must, notes)
        # Every file the notes name exists in the release workflow's list.
        names = set(release_file_names("1.2.3"))
        found = [f for f in re.findall(r"`((?:COAL|coal)[\w.-]*\.(?:exe|zip|whl|tar\.gz))`", notes) if "..." not in f]
        self.assertGreaterEqual(len(found), 6)
        for f in found:
            self.assertIn(f, names)

    def test_changelog_section(self):
        rn = _load("release_notes")
        with tempfile.TemporaryDirectory() as t:
            p = os.path.join(t, "CHANGELOG.md")
            with open(p, "w", encoding="utf-8") as fh:
                fh.write("# Changelog\n\n## [1.2.3] - 2026-10-18\n\n- New thing.\n\n## v1.2.2\n\n- Old thing.\n")
            self.assertEqual(rn.changelog_section("1.2.3", p), "- New thing.")
            self.assertEqual(rn.changelog_section("1.2.2", p), "- Old thing.")
            self.assertIsNone(rn.changelog_section("1.2", p))
            self.assertIsNone(rn.changelog_section("1.2.3", os.path.join(t, "missing.md")))


class TestReleaseWorkflow(unittest.TestCase):
    def test_parses_and_triggers(self):
        wf = _workflow("release.yml")
        on = wf.get("on", wf.get(True))  # YAML 1.1 reads a bare `on` key as True
        self.assertEqual(on["push"]["tags"], ["v*"])
        self.assertIn("workflow_dispatch", on)
        self.assertEqual(set(wf["jobs"]), {"python", "windows", "release"})

    def test_permissions_minimal(self):
        wf = _workflow("release.yml")
        self.assertEqual(wf["permissions"], {"contents": "read"})
        for name, job in wf["jobs"].items():
            if name == "release":
                self.assertEqual(job["permissions"], {"contents": "write"})
            else:
                self.assertNotIn("permissions", job, name)

    def test_actions_pinned_to_commits(self):
        for fname in ("release.yml",) + (("pages.yml",) if os.path.exists(os.path.join(WORKFLOWS, "pages.yml"))
                                          else ()):
            text = _read(".github", "workflows", fname)
            uses = re.findall(r"uses:\s*(\S+)", text)
            self.assertTrue(uses, fname)
            for u in uses:
                self.assertRegex(u, r"^[\w.-]+/[\w./-]+@[0-9a-f]{40}$", f"{fname}: {u} is not pinned to a commit")

    def test_windows_job_builds_and_checks(self):
        wf = _workflow("release.yml")
        runs = "\n".join(s.get("run", "") for s in wf["jobs"]["windows"]["steps"])
        self.assertEqual(wf["jobs"]["windows"]["runs-on"], "windows-latest")
        self.assertIn(r"desktop\build.ps1", runs)
        self.assertRegex(runs, r"choco install innosetup --version=\d+\.\d+\.\d+")
        self.assertIn('"-m", "coal", "--version"', runs)
        self.assertIn("Get-FileHash", runs)  # frozen vs Python output comparison
        for name in ("COAL-Setup-$v.exe", "COAL-Setup-windows-x64.exe", "COAL-$v-windows-x64.zip",
                     "COAL-windows-x64.zip"):
            self.assertIn(name, runs)

    def test_every_listed_file_is_produced(self):
        """Each name in the release check list is created by some build step (or is the checksum file)."""
        wf = _workflow("release.yml")
        made = "\n".join(s.get("run", "") for j in ("python", "windows") for s in wf["jobs"][j]["steps"])
        made = made.replace("$v", "${V}").replace("${env:VERSION}", "${V}")
        for name in release_file_names("${V}"):
            if name == "SHA256SUMS.txt":
                continue
            self.assertTrue(name in made or name in ("coal-${V}-py3-none-any.whl", "coal-${V}.tar.gz")
                            and "python -m build" in made, f"no build step makes {name}")
        sums = _steps(wf, "release")["SHA-256 of every file"]["run"]
        self.assertIn("> SHA256SUMS.txt", sums)


if __name__ == "__main__":
    unittest.main()
