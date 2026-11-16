"""The download website (site/), the release/pages workflows behind it, and the scripts that adapt the deployed copy to the latest release (or to none yet)."""
import importlib.util
import os
import re
import shutil
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


def platform_wheel_tags():
    """The platform tags of the engine wheels with coal-cpp inside, from the cpp-wheels job's matrix."""
    return [m["tag"] for m in _workflow("release.yml")["jobs"]["cpp-wheels"]["strategy"]["matrix"]["include"]]


def release_file_names(version="9.9.9"):
    """The exact file list release.yml checks (the files job, step "Check the release files"), with ${V} filled in
    and, when the version has an R form (not a prerelease), coalr's files with ${R_V} = r_version(version)."""
    run = _steps(_workflow("release.yml"), "files")["Check the release files"]["run"]
    block = re.search(r"expected=\((.*?)\)", run, re.S).group(1)
    names = re.findall(r'"([^"]+)"', block)
    r_v = _load("release_version").r_version(version.replace("${V}", "9.9.9"))
    if r_v:
        extra = re.search(r'if \[ -n "\$R_V" \]; then expected\+=\((.*?)\); fi', run, re.S)
        names += re.findall(r'"([^"]+)"', extra.group(1))
    r_name = "${R_V}" if version == "${V}" else (r_v or "")
    return [n.replace("${V}", version).replace("${R_V}", r_name) for n in names]


class TestReleaseScripts(unittest.TestCase):
    def test_version_check(self):
        rv = _load("release_version")
        init, proj = rv.source_versions()
        self.assertEqual(init, proj, "src/coal/__init__.py and pyproject.toml disagree")
        # one version for both packages: the app's pyproject and the coal-align version it pins
        app, pin = rv.app_versions()
        self.assertEqual(app, init, "apps/coal-app/pyproject.toml has another version than the engine")
        self.assertEqual(pin, init, "apps/coal-app/pyproject.toml must depend on coal-align==<the same version>")
        self.assertEqual(rv.cpp_versions(), (init, init), "cpp/src/main.cpp and report.cpp VERSION must be the engine's")
        if rv.DEV_RE.match(init):   # a development tree (0.6.0.dev0) has no tag: the release commit sets 0.6.0 first
            for tag in ("v" + init, "v" + rv.base_version(init)):
                with self.assertRaisesRegex(ValueError, "development version|is not v<major>", msg=tag):
                    rv.check(tag)
        else:
            self.assertEqual(rv.check("v" + init), (init, bool(re.search(r"(a|b|rc)\d+$", init))))
        self.assertEqual([rv.base_version(v) for v in ("0.6.0.dev0", "0.5.0rc1", "0.5.0")], ["0.6.0", "0.5.0", "0.5.0"])
        self.assertEqual([rv.r_version(v) for v in ("0.6.0.dev0", "0.6.0.dev3", "0.6.0", "0.6.1rc1", "1.0.0b2")],
                         ["0.6.0.9000", "0.6.0.9003", "0.6.0", None, None])
        self.assertEqual(rv.description_version(), rv.r_version(init) or rv.description_version())
        for bad in ("v9.9.9", init, "v" + init + "x", "", "release-" + init):
            with self.assertRaises(ValueError, msg=bad):
                rv.check(bad)
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "release_version.py"), "v9.9.9"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 1)
        self.assertIn("bump all six", r.stderr)
        self.assertIn("cpp/src/main.cpp", r.stderr)
        self.assertIn("apps/coal-app/pyproject.toml", r.stderr)

    def test_dry_run_and_release_conditions(self):
        """--dry-run (release.yml dispatched with publish off on a branch): the source's version, once the six
        places agree; a tag is also refused while CHANGELOG.md has "Unreleased" entries or coalr's DESCRIPTION has
        another version (tested on a copy of the files the check reads)."""
        rv = _load("release_version")
        init = rv.source_versions()[0]
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "out.txt")
            r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "release_version.py"), "--dry-run",
                                "main", "--github-output"], capture_output=True, text=True, timeout=60,
                               env={**os.environ, "GITHUB_OUTPUT": out})
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(r.stdout, init + "\n")
            self.assertIn(f"dry run of main: version {init}", r.stderr)
            with open(out, encoding="utf-8") as fh:
                self.assertEqual(fh.read(), f"version={init}\nprerelease=true\nr_version={rv.r_version(init) or ''}\n")
            # a copy of the tree's version files at release 1.2.3, then each thing the check refuses
            for rel in ("src/coal/__init__.py", "pyproject.toml", "apps/coal-app/pyproject.toml", "cpp/src/main.cpp",
                        "cpp/src/report.cpp", "CHANGELOG.md", "packaging/r/coalr/DESCRIPTION",
                        "packaging/r/coalr/NEWS.md"):
                dst = os.path.join(t, "tree", *rel.split("/"))
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                with open(os.path.join(ROOT, *rel.split("/")), encoding="utf-8") as fh:
                    text = fh.read().replace(init, "1.2.3")
                if rel == "CHANGELOG.md":
                    text = "# Changelog\n\n## Unreleased\n\n## 1.2.3\n\n- x\n"
                if rel.endswith("DESCRIPTION"):
                    text = re.sub(r"(?m)^Version: .*$", "Version: 1.2.3", text)
                if rel.endswith("NEWS.md"):
                    text = "# coalr 1.2.3\n\n- x\n\n# coalr 1.2.2\n"
                with open(dst, "w", encoding="utf-8") as fh:
                    fh.write(text)
            tree = os.path.join(t, "tree")
            self.assertEqual(rv.check("v1.2.3", tree), ("1.2.3", False))
            self.assertEqual(rv.dry_run(tree), "1.2.3")
            log = os.path.join(tree, "CHANGELOG.md")
            with open(log, "w", encoding="utf-8") as fh:
                fh.write("# Changelog\n\n## Unreleased\n\n### Fixed\n\n- y\n\n## 1.2.3\n\n- x\n")
            with self.assertRaisesRegex(ValueError, "entries under \"## Unreleased\""):
                rv.check("v1.2.3", tree)
            self.assertEqual(rv.dry_run(tree), "1.2.3")          # a dry run publishes nothing: not refused
            with open(log, "w", encoding="utf-8") as fh:
                fh.write("# Changelog\n\n## 1.2.3\n\n- x\n")
            desc = os.path.join(tree, "packaging", "r", "coalr", "DESCRIPTION")
            with open(desc, encoding="utf-8") as fh:
                text = fh.read()
            with open(desc, "w", encoding="utf-8") as fh:
                fh.write(text.replace("Version: 1.2.3", "Version: 1.2.2.9000"))
            with self.assertRaisesRegex(ValueError, "DESCRIPTION says Version: 1.2.2.9000, not 1.2.3"):
                rv.check("v1.2.3", tree)
            with open(desc, "w", encoding="utf-8") as fh:
                fh.write(text)
            news = os.path.join(tree, "packaging", "r", "coalr", "NEWS.md")
            with open(news, "w", encoding="utf-8") as fh:
                fh.write("# coalr (development version)\n\n- y\n\n# coalr 1.2.2\n")
            with self.assertRaisesRegex(ValueError, 'NEWS.md starts with "# coalr \\(development version\\)", not "# coalr 1.2.3"'):
                rv.check("v1.2.3", tree)
            with open(news, "w", encoding="utf-8") as fh:
                fh.write("# coalr 1.2.3\n")
            self.assertEqual(rv.check("v1.2.3", tree), ("1.2.3", False))
            with open(desc, "w", encoding="utf-8") as fh:
                fh.write(text.replace("Version: 1.2.3", "Version: 1.2.2.9000"))
            with open(os.path.join(tree, "src", "coal", "__init__.py"), "a", encoding="utf-8") as fh:
                fh.write('__version__ = "1.2.4"\n')   # the last assignment is not read: the first one is
            with open(os.path.join(tree, "pyproject.toml"), encoding="utf-8") as fh:
                proj = fh.read()
            with open(os.path.join(tree, "pyproject.toml"), "w", encoding="utf-8") as fh:
                fh.write(proj.replace('version = "1.2.3"', 'version = "1.2.4"', 1))
            with self.assertRaisesRegex(ValueError, "dry run: the versions disagree"):
                rv.dry_run(tree)

    def test_prerelease_tags(self):
        rv = _load("release_version")
        self.assertTrue(rv.TAG_RE.match("v1.0.0rc1"))
        self.assertFalse(rv.TAG_RE.match("v1.0.0-rc1"))

    def test_notes(self):
        rn = _load("release_notes")
        notes = rn.fixed_part("1.2.3")
        for must in ("not code-signed", "Windows protected your PC", "Smart App Control", "Run anyway",
                     "pip install", "coal_align-1.2.3-py3-none-any.whl", "releases/download/v1.2.3/", "SHA256SUMS.txt",
                     "Get-FileHash", "sha256sum"):
            self.assertIn(must, notes)
        # Every file the notes name exists in the release workflow's list.
        names = set(release_file_names("1.2.3"))
        found = [f for f in re.findall(r"`((?:COAL|coal)[\w.-]*\.(?:exe|zip|whl|tar\.gz))`", notes) if "..." not in f]
        self.assertGreaterEqual(len(found), 6)
        for f in found:
            self.assertIn(f, names)
        # and every file of the list is in the notes' Downloads table, coalr's too (none for a prerelease)
        for version in ("1.2.3", "1.2.3rc1"):
            notes = rn.fixed_part(version)
            table = notes[notes.index("## Downloads"):notes.index("The version-free names")]
            names = release_file_names(version)
            for f in names:
                self.assertIn(f"`{f}`", table, f"{version}: the release notes do not list {f}")
            self.assertEqual("coalr_" in table, version == "1.2.3", version)
        self.assertIn("releases/download/v1.2.3/coalr_1.2.3.tar.gz", rn.fixed_part("1.2.3"))

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

    def test_notes_fit_a_github_release(self):
        """GitHub refuses a release body over 125,000 characters, which 0.6.0's whole CHANGELOG section is: the
        notes then carry its summary and a link to CHANGELOG.md at the tag."""
        rn = _load("release_notes")
        s = rn.summary_of("**In short.** Big.\n\n- one\n\n### Added\n\n- x\n", "1.2.3")
        self.assertTrue(s.startswith("**In short.** Big.\n\n- one\n\nThe full list"), s)
        self.assertNotIn("### Added", s)
        self.assertIn("https://github.com/plumat-coal/COAL/blob/v1.2.3/CHANGELOG.md", s)
        v = _load("release_version").base_version(_load("release_version").source_versions()[0])
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "notes.md")
            rc = rn.main([v, out])
            text = open(out, encoding="utf-8").read()
        self.assertIn(rc, (0, rn.NO_CHANGELOG))
        self.assertLess(len(text), 125_000)
        if rc == 0:
            self.assertIn(f"## Changes in {v}", text)


class TestReleaseWorkflow(unittest.TestCase):
    def test_parses_and_triggers(self):
        wf = _workflow("release.yml")
        on = wf.get("on", wf.get(True))  # YAML 1.1 reads a bare `on` key as True
        self.assertEqual(on["push"]["tags"], ["v*"])
        self.assertIn("workflow_dispatch", on)
        self.assertEqual(set(wf["jobs"]), {"python", "cpp-wheels", "windows", "coalr", "files", "release", "pypi"})

    def test_coalr_job(self):
        """coalr is a release asset (REVIEW-COAL-3 H1): built from packaging/r/coalr after the vendored-source and
        version checks, installed and run, uploaded versioned and version-free; skipped for a prerelease, whose
        release then has no coalr files; the files and release jobs still run then."""
        wf = _workflow("release.yml")
        job = wf["jobs"]["coalr"]
        self.assertEqual(job["needs"], "python")
        self.assertEqual(job["if"], "${{ needs.python.outputs.r_version != '' }}")
        self.assertEqual(job["env"]["R_V"], "${{ needs.python.outputs.r_version }}")
        self.assertEqual(wf["jobs"]["python"]["outputs"]["r_version"], "${{ steps.v.outputs.r_version }}")
        runs = "\n".join(s.get("run", "") for s in job["steps"])
        for must in ("python scripts/sync_r_sources.py --check", 'grep -qx "Version: ${R_V}" packaging/r/coalr/DESCRIPTION',
                     "R CMD build --no-manual packaging/r/coalr", 'R CMD INSTALL "coalr_${R_V}.tar.gz"',
                     'packageVersion("coalr")', "coal_dnds_pairs(", 'cp "coalr_${R_V}.tar.gz" release/coalr_src.tar.gz'):
            self.assertIn(must, runs)
        uses = [s["uses"].split("@")[0] for s in job["steps"] if "uses" in s]
        self.assertEqual(uses, ["actions/checkout", "actions/setup-python", "r-lib/actions/setup-r",
                                "r-lib/actions/setup-pandoc", "r-lib/actions/setup-r-dependencies",
                                "actions/upload-artifact"])
        self.assertEqual(job["steps"][-1]["with"]["name"], "dist-coalr")
        # the same r-lib actions (SHAs) as r-package.yml
        rpkg = {u.split("@")[0]: u for u in re.findall(r"uses: (\S+)", _read(".github", "workflows", "r-package.yml"))}
        for st in job["steps"]:
            if st.get("uses", "").startswith("r-lib/"):
                self.assertEqual(st["uses"], rpkg[st["uses"].split("@")[0]])
        self.assertIn("coalr", wf["jobs"]["files"]["needs"])
        self.assertEqual(wf["jobs"]["release"]["needs"], ["python", "files"])
        for j in ("files", "release", "pypi"):
            self.assertIn("!cancelled() && !failure()", wf["jobs"][j]["if"], j)   # a skipped coalr does not skip them
        self.assertIn("coalr_1.2.3.tar.gz", release_file_names("1.2.3"))
        self.assertIn("coalr_src.tar.gz", release_file_names("1.2.3"))
        self.assertIn("coalr_1.2.3.9000.tar.gz", release_file_names("1.2.3.dev0"))
        self.assertFalse([n for n in release_file_names("1.2.3rc1") if n.startswith("coalr")])
        check = _steps(wf, "files")["Check the release files"]
        self.assertEqual(check["env"]["PRE"], "${{ needs.python.outputs.prerelease }}")
        self.assertIn('if [ "$PRE" = false ] && [ -z "$R_V" ]; then', check["run"])   # a release needs coalr

    def test_bash_steps_parse(self):
        """Every bash step of release.yml and release-smoke.yml is valid shell (bash -n): a quote inside the R
        snippet of the coalr job ended its single-quoted script on the first dry run."""
        bash = shutil.which("bash")
        if not bash:
            self.skipTest("no bash")
        for fname in ("release.yml", "release-smoke.yml"):
            for name, job in _workflow(fname)["jobs"].items():
                default = job.get("defaults", {}).get("run", {}).get("shell")
                for st in job["steps"]:
                    if "run" in st and (st.get("shell") or default or "bash") == "bash":
                        r = subprocess.run([bash, "-n"], input=st["run"], text=True, capture_output=True, timeout=60)
                        self.assertEqual(r.returncode, 0, f"{fname} {name} / {st.get('name')}: {r.stderr}")

    def test_dry_run_and_existing_release(self):
        """The dry-run version check is scripts/release_version.py --dry-run (tested above), not inline Python; an
        existing release gets the same prerelease flag as a new one."""
        wf = _workflow("release.yml")
        for j in ("python", "windows"):
            run = _steps(wf, j)["Tag matches the version in the source"]["run"]
            self.assertIn('python scripts/release_version.py --dry-run "$TAG" --github-output', run, j)
            self.assertIn('python scripts/release_version.py "$TAG" --github-output', run, j)
            self.assertNotIn("<<", run, j)
        run = _steps(wf, "release")["Create or update the GitHub Release"]["run"]
        self.assertIn('gh release edit "$TAG" --title "COAL $V" --notes-file notes.md "--prerelease=${PRE}"', run)
        self.assertIn('if [ "$PRE" = "true" ]; then extra+=(--prerelease); fi', run)
        # the release publishes exactly what the files job checked
        downloads = [st["with"] for st in wf["jobs"]["release"]["steps"]
                     if st.get("uses", "").startswith("actions/download-artifact@")]
        self.assertEqual(downloads, [{"name": "release-files", "path": "release"}])

    def test_pypi_job_is_opt_in_and_uses_trusted_publishing(self):
        job = _workflow("release.yml")["jobs"]["pypi"]
        self.assertTrue(job["if"].endswith("&& vars.PYPI_PUBLISH == 'true' }}"), job["if"])  # off until the variable is set
        self.assertEqual(job["environment"], "pypi")
        self.assertEqual(job["permissions"], {"id-token": "write", "contents": "read"})  # no API token
        self.assertIn("release", job["needs"])  # the GitHub Release exists before PyPI gets the files
        uses = [st["uses"] for st in job["steps"] if "uses" in st]
        self.assertTrue(any(u.startswith("pypa/gh-action-pypi-publish@") and "#" not in u and len(u.split("@")[1]) == 40
                            for u in uses), uses)  # pinned to a commit

    def test_permissions_minimal(self):
        wf = _workflow("release.yml")
        self.assertEqual(wf["permissions"], {"contents": "read"})
        for name, job in wf["jobs"].items():
            if name == "release":
                self.assertEqual(job["permissions"], {"contents": "write"})
            elif name == "pypi":
                self.assertEqual(job["permissions"], {"id-token": "write", "contents": "read"})
            else:
                self.assertNotIn("permissions", job, name)

    def test_ci_caches_pip_and_browsers(self):
        for fname in ("ci.yml", "nightly.yml"):
            for name, job in _workflow(fname)["jobs"].items():
                for st in job["steps"]:
                    if st.get("uses", "").startswith("actions/setup-python@"):
                        self.assertEqual(st["with"].get("cache"), "pip", f"{fname}: {name}")
                        self.assertIn("pyproject.toml", st["with"].get("cache-dependency-path", ""), f"{fname}: {name}")
        wf = _workflow("ci.yml")
        cache = [st for st in wf["jobs"]["explainer"]["steps"] if st.get("uses", "").startswith("actions/cache@")]
        self.assertEqual(len(cache), 1)
        self.assertEqual(cache[0]["with"]["path"], wf["jobs"]["explainer"]["env"]["PLAYWRIGHT_BROWSERS_PATH"])
        self.assertIn("playwright-1.48.0", cache[0]["with"]["key"])  # same pin as the install step

    def test_actions_pinned_to_commits(self):
        # ci.yml too (it uses actions/cache), and the nightly fuzz workflow
        for fname in ("release.yml", "ci.yml", "nightly.yml", "r-package.yml") + (("pages.yml",) if os.path.exists(os.path.join(WORKFLOWS, "pages.yml"))
                                          else ()):
            text = _read(".github", "workflows", fname)
            uses = re.findall(r"uses:\s*(\S+)", text)
            self.assertTrue(uses, fname)
            for u in uses:
                self.assertRegex(u, r"^[\w.-]+/[\w./-]+@[0-9a-f]{40}$", f"{fname}: {u} is not pinned to a commit")

    def test_nightly_workflow_is_scheduled_read_only_and_bounded(self):
        wf = _workflow("nightly.yml")
        on = wf.get("on", wf.get(True))  # PyYAML reads the bare key `on` as True
        self.assertEqual(on["schedule"], [{"cron": "17 3 * * *"}])
        self.assertIn("seed", on["workflow_dispatch"]["inputs"])
        self.assertFalse(on["workflow_dispatch"]["inputs"]["seed"]["required"])
        self.assertEqual(wf["permissions"], {"contents": "read"})
        self.assertIn("group", wf["concurrency"])
        self.assertEqual(wf["env"]["FUZZ_SEED"], "${{ inputs.seed || github.run_number }}")
        self.assertEqual(set(wf["jobs"]), {"pal2nal-fuzz", "engine-fuzz", "kaks-properties", "sanitizers", "tsan", "cpp-extras",
                                         "msa-tools", "tutorials"})
        for name, job in wf["jobs"].items():
            self.assertIsInstance(job.get("timeout-minutes"), int, name)
            self.assertNotIn("permissions", job, name)
            uploads = [st for st in job["steps"] if st.get("uses", "").startswith("actions/upload-artifact@")]
            self.assertEqual(len(uploads), 1, name)
            self.assertEqual(uploads[0]["if"], "failure()", name)  # the kept cases, only when something differed
            self.assertIn("retention-days", uploads[0]["with"], name)

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
        made = "\n".join(s.get("run", "") for j in ("python", "windows", "coalr") for s in wf["jobs"][j]["steps"])
        made = made.replace("$v", "${V}").replace("${env:VERSION}", "${V}")
        # the cpp-wheels job names its wheel through the matrix tag
        wheels = "\n".join(s.get("run", "") for s in wf["jobs"]["cpp-wheels"]["steps"])
        for tag in platform_wheel_tags():
            made += "\n" + wheels.replace("${TAG_NAME}", tag)
        for name in release_file_names("${V}"):
            if name == "SHA256SUMS.txt":
                continue
            self.assertTrue(name in made or name in ("coal_align-${V}-py3-none-any.whl", "coal_align-${V}.tar.gz")
                            and "python -m build" in made, f"no build step makes {name}")
        sums = _steps(wf, "files")["SHA-256 of every file"]["run"]
        self.assertIn("> SHA256SUMS.txt", sums)

    def test_platform_wheels_job(self):
        """One engine wheel with coal-cpp per platform (issue #16): built from the pure wheel by
        packaging/wheels/build_platform_wheel.py, checked (auditwheel on Linux, both slices on macOS, an install
        whose `coal --cpp-path` is the bundled binary and whose output equals the Python engine's) and uploaded for
        the release and PyPI jobs; the pure wheel stays in the release as the fallback."""
        wf = _workflow("release.yml")
        job = wf["jobs"]["cpp-wheels"]
        self.assertEqual(job["needs"], "python")  # the pure wheel comes from there
        tags = platform_wheel_tags()
        self.assertEqual(tags, ["manylinux_2_28_x86_64", "manylinux_2_28_aarch64", "macosx_11_0_universal2",
                                "win_amd64"])
        rn = _load("release_notes")
        self.assertEqual(list(rn.PLATFORM_TAGS), tags)  # the notes list the same wheels
        by_tag = {m["tag"]: m for m in job["strategy"]["matrix"]["include"]}
        for tag in ("manylinux_2_28_x86_64", "manylinux_2_28_aarch64"):
            self.assertEqual(by_tag[tag]["image"], f"quay.io/pypa/{tag}", tag)
        self.assertEqual(by_tag["manylinux_2_28_aarch64"]["os"], "ubuntu-24.04-arm")  # native, no QEMU
        self.assertTrue(by_tag["macosx_11_0_universal2"]["os"].startswith("macos"))
        self.assertTrue(by_tag["win_amd64"]["os"].startswith("windows"))
        runs = "\n".join(s.get("run", "") for s in job["steps"])
        self.assertIn("packaging/wheels/build_linux_binary.sh", runs)
        self.assertIn("packaging/wheels/build_platform_wheel.py", runs)
        self.assertIn("--audit", runs)  # auditwheel confirms the manylinux tag
        self.assertIn('-DCMAKE_OSX_ARCHITECTURES="arm64;x86_64"', runs)
        self.assertIn("lipo -archs", runs)
        self.assertIn("-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded", runs)  # no VC redistributable needed
        self.assertIn("-DCOAL_WARNINGS_AS_ERRORS=ON", runs)
        self.assertIn("--cpp-path", runs)
        self.assertIn("COAL_ENGINE=cpp", runs)
        self.assertIn("cmp ", runs)
        uploads = [st for st in job["steps"] if st.get("uses", "").startswith("actions/upload-artifact@")]
        self.assertEqual(len(uploads), 1)
        self.assertEqual(uploads[0]["with"]["name"], "dist-wheel-${{ matrix.tag }}")
        for j in ("files", "pypi"):
            self.assertIn("cpp-wheels", wf["jobs"][j]["needs"], j)
        names = release_file_names("1.2.3")
        self.assertIn("coal_align-1.2.3-py3-none-any.whl", names)  # the fallback wheel is still released
        for tag in tags:
            self.assertIn(f"coal_align-1.2.3-py3-none-{tag}.whl", names)
        pypi = "\n".join(s.get("run", "") for s in wf["jobs"]["pypi"]["steps"])
        self.assertIn("coal_align-${V}-py3-none-*.whl", pypi)
        self.assertIn("-eq 5", pypi)  # the pure wheel and the four platform wheels go to PyPI


class TestPagesWorkflow(unittest.TestCase):
    def test_deploys_site_on_main_with_minimal_permissions(self):
        wf = _workflow("pages.yml")
        on = wf.get("on", wf.get(True))
        self.assertEqual(on["push"]["branches"], ["main"])
        self.assertIn("site/**", on["push"]["paths"])
        self.assertIn("workflow_dispatch", on)
        self.assertEqual(wf["permissions"], {"contents": "read"})
        self.assertNotIn("permissions", wf["jobs"]["build"])
        self.assertEqual(wf["jobs"]["deploy"]["permissions"], {"pages": "write", "id-token": "write"})
        self.assertEqual(wf["jobs"]["deploy"]["environment"]["name"], "github-pages")
        uses = [s.get("uses", "") for j in wf["jobs"].values() for s in j["steps"]]
        self.assertTrue(any(u.startswith("actions/upload-pages-artifact@") for u in uses))
        self.assertTrue(any(u.startswith("actions/deploy-pages@") for u in uses))
        build = "\n".join(s.get("run", "") for s in wf["jobs"]["build"]["steps"])
        self.assertIn("cp -R site/. _site/", build)
        self.assertIn("tests/test_docs_site.py", build)


# Several pages, a self-hosted font and its licence, a WebP screenshot, logo, mark and icons.
PAGES = ("index.html", "install.html", "guide.html", "differences.html", "cite.html")
FONTS = tuple(f"assets/fonts/Literata-{w}.woff2" for w in ("Regular", "Semibold", "Italic"))
SITE_FILES = PAGES + FONTS + ("assets/fonts/OFL.txt",) + ("index.html", "style.css", "main.js", "sw.js", "favicon.svg", "explainer.html", "img/broken-gene.webp",
              "assets/coal-mark.svg", "assets/coal-mark-dark.svg", "assets/coal-logo.png", "assets/coal-logo@2x.png", "assets/coal-logo-dark.png",
              "assets/coal-logo-dark@2x.png", "assets/favicon-32.png", "assets/apple-touch-icon.png",
              "assets/icon-512.png", "assets/coal-train.svg")
TEXT_EXT = (".html", ".css", ".js", ".svg", ".txt", ".json", ".webmanifest")
# Not requests: the SVG namespace, the local app's address shown as text for people to open, and the DOI links of
# the cited papers on cite.html (links a reader may follow; the page loads nothing from them. They are the DOIs of
# other people's papers: the rule "no DOI for COAL itself" below is unaffected).
ALLOWED_NON_REPO_URLS = ("http://www.w3.org/2000/svg", "http://127.0.0.1:7860", "https://doi.org/10.")


# The font's licence names its project and licence sites; no page loads or links them.
NOT_LOADED = (os.path.join("assets", "fonts", "OFL.txt"),)


def site_text_files():
    out = []
    for d, _, files in os.walk(SITE):
        for f in files:
            if os.path.relpath(os.path.join(d, f), SITE) in NOT_LOADED:
                continue
            if f.endswith(TEXT_EXT):
                out.append(os.path.join(d, f))
    return sorted(out)


def site_links(html):
    return re.findall(r'(?:href|src)\s*=\s*"([^"]*)"', html)


class TestSiteFiles(unittest.TestCase):
    # content is spread over index, install, guide, differences and cite pages. Checks that need
    # a string somewhere on the site look at all pages (self.all); per-page checks loop over self.pages.
    def setUp(self):
        self.index = _read("site", "index.html")
        self.pages = {p: _read("site", p) for p in PAGES}
        self.all = "\n".join(self.pages.values())

    def test_files_exist(self):
        for f in SITE_FILES + ("README.md",):
            self.assertTrue(os.path.isfile(os.path.join(SITE, f)), f"site/{f} is missing")

    def test_no_external_urls(self):
        """Only links to github.com/plumat-coal/COAL; nothing else may be referenced (no CDNs, fonts, trackers)."""
        for path in site_text_files():
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
            rel = os.path.relpath(path, ROOT)
            for url in re.findall(r"""(?:https?:)?//[^\s"'<>()]+""", text):
                url = url.rstrip(".,;:")  # sentence punctuation after a URL in running text
                if url in ALLOWED_NON_REPO_URLS or url.startswith(ALLOWED_NON_REPO_URLS):
                    continue
                if url.startswith("//") and not url.startswith("//www") and "." not in url.split("/")[2]:
                    continue  # a JS comment ("// ..."), not a protocol-relative URL
                self.assertRegex(url, r"^https://github\.com/plumat-coal/COAL(?:[/#?]|$)", f"{rel}: external URL {url}")
            for link in site_links(text):
                self.assertFalse(link.startswith("//"), f"{rel}: protocol-relative link {link}")
                if not link.startswith(("https://", "#", "mailto:")):
                    self.assertFalse(link.startswith("/"), f"{rel}: {link} must be relative (site may move to a custom domain)")

    def test_internal_links_resolve(self):
        # every page, and page#anchor links must find the anchor on that page
        ids = {p: set(re.findall(r'\sid="([^"]+)"', h)) for p, h in self.pages.items()}
        for page, html in self.pages.items():
            for link in site_links(html):
                if link.startswith("#"):
                    self.assertIn(link[1:], ids[page], f"{page}: anchor {link} has no target")
                elif not link.startswith("https://"):
                    target = link.split("#")[0].split("?")[0] or page
                    if target in ("./", "."):
                        target = "index.html"
                    self.assertTrue(os.path.isfile(os.path.join(SITE, target)), f"{page}: link {link}: site/{target} missing")
                    if "#" in link and target in ids:
                        self.assertIn(link.split("#")[1], ids[target], f"{page}: link {link} has no target")
        sw = _read("site", "sw.js")
        for f in re.search(r"FILES\s*=\s*\[(.*?)\]", sw, re.S).group(1).replace('"', "").split(","):
            f = f.strip()
            if f not in ("./", ""):
                self.assertTrue(os.path.isfile(os.path.join(SITE, f)), f"sw.js caches missing file {f}")

    def test_download_links_match_release_files(self):
        names = set(release_file_names())
        # all pages
        latest = re.findall(r"https://github\.com/plumat-coal/COAL/releases/latest/download/([^\"'\s<]+)", self.all)
        self.assertGreaterEqual(len(set(latest)), 4, "expected installer, zip, sdist and checksum links")
        for f in latest:
            self.assertIn(f, names, f"site links to {f}, which the release workflow does not upload")
            self.assertNotRegex(f, r"\d+\.\d+", f"{f}: latest/download links need version-free names")
        self.assertNotIn("/releases/download/", self.all, "versioned download links go stale; use releases/latest")
        self.assertIn('id="win-download" href="https://github.com/plumat-coal/COAL/releases/latest/download/'
                      'COAL-Setup-windows-x64.exe"', self.index)

    def test_pip_commands_use_the_cache_safe_sdist(self):
        cmds = re.findall(r'"coal-align(?:\[gui\])? @ ([^"]+)"', self.all)  # all pages
        self.assertGreaterEqual(len(cmds), 4)
        for url in cmds:
            self.assertTrue(url.endswith("/releases/latest/download/coal_src.tar.gz"), url)
        # the app (the second package) comes from its own version-free sdist, always together with the engine's
        app = re.findall(r'"coal-align-app(?:\[gui\])? @ ([^"]+)"', self.all)
        self.assertGreaterEqual(len(app), 4)
        for url in app:
            self.assertTrue(url.endswith("/releases/latest/download/coal_app_src.tar.gz"), url)
        for line in re.findall(r"[^\n]*coal-align-app\[gui\] @[^\n]*", self.all):
            self.assertIn('"coal-align @ ', line, f"the app installed without the engine's sdist: {line}")
        # pip caches wheels built from URLs whose file name looks like name-version; the aliases must not.
        self.assertNotIn("-", "coal_src")
        self.assertNotIn("-", "coal_app_src")

    def test_claims_are_in_the_readme(self):
        """Every measured number on the page (data-check) is quoted from README.md, docs/ or tests/README.md, which
        test_docs_numbers guards (the test counts are in tests/README.md)."""
        docs = [os.path.join("docs", f) for f in sorted(os.listdir(os.path.join(ROOT, "docs"))) if f.endswith(".md")]
        readme = re.sub(r"\s+", " ", " ".join(_read(*p.split(os.sep)) for p in ["README.md", *docs, "tests/README.md"]))
        readme = readme.replace("**", "")
        claims = re.findall(r'data-check="([^"]+)">([^<]*)<', self.all)  # all pages
        self.assertGreaterEqual(len(claims), 10)
        for value, shown in claims:
            self.assertEqual(value, shown, "data-check value and the text shown differ")
            self.assertIn(value, readme, f"site claims {value!r}, which README.md does not say")
        self.assertIn("Measured on", self.pages["differences.html"])

    def test_version_is_current(self):
        """The page names the current version; it must match the package."""
        init = re.search(r'__version__ = "([^"]+)"', _read("src", "coal", "__init__.py")).group(1)
        init = _load("release_version").base_version(init)   # a development tree shows the release it becomes
        shown = re.search(r'<span id="version" data-version="([^"]+)">([^<]+)</span>', self.index)
        self.assertIsNotNone(shown, "version span missing")
        self.assertEqual(shown.group(1), shown.group(2))
        self.assertEqual(shown.group(1), init, "site/index.html names another version than src/coal/__init__.py")
        # the citation (on cite.html) names the version too
        ref = re.findall(r"<span data-version-ref>([^<]+)</span>", self.all)
        self.assertEqual(ref, [init], "site/cite.html: the citation's data-version-ref span names another version")

    def test_pal2nal_relationship_and_citation(self):
        """The qualified PAL2NAL wording, the licences and the reference."""
        # the documented-exceptions wording (an "on well-formed input" claim would be false)
        for must in ("COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate "
                     "differences", "GPL-2.0", "contains none of PAL2NAL's code",
                     "Suyama M, Torrents D, Bork P (2006)", "34:W609"):
            self.assertIn(must, self.all)  # differences.html and cite.html

    def test_size_budget_with_images(self):
        """The whole site, images included, stays under 400 KB (README.md excluded: not deployed)."""
        total = sum(os.path.getsize(os.path.join(d, f)) for d, _, files in os.walk(SITE) for f in files
                    if f != "README.md")
        # the results screenshot is sharp on high-resolution screens (1600 px wide, ~175 KB), so 400 KB
        self.assertLess(total, 400 * 1024, f"site is {total} bytes (budget 400 KB)")

    def test_cite_table_wraps_and_stacks_on_a_phone(self):
        """The "what to cite" table has prose in both columns: its own class, not the file table's nowrap first
        column (on a 390 px screen only the first column was visible), and stacked cells below 30rem."""
        cite = self.pages["cite.html"]
        self.assertEqual(re.findall(r'<table class="([^"]+)"', cite), ["files refs"])
        for page, html in self.pages.items():
            if page != "cite.html":
                self.assertNotIn("refs", re.findall(r'<table class="([^"]+)"', html), page)
        css = re.sub(r"\s+", " ", _read("site", "style.css"))
        self.assertIn(".refs td:first-child { font: inherit; white-space: normal; }", css)
        self.assertLess(css.index(".files td:first-child {"), css.index(".refs td:first-child {"))   # the override wins
        phone = re.search(r"@media \(max-width: 30rem\) \{ \.refs \{ min-width: 0; \}.*?\} \}", css)
        self.assertIsNotNone(phone, "no phone rule for .refs")
        self.assertIn(".refs td { display: block;", phone.group(0))
        self.assertNotIn("display: none", phone.group(0))   # the header row stays for screen readers

    def test_unsigned_and_smart_app_control_explained(self):
        for must in ("Windows protected your PC", "More info", "Run anyway", "Unknown publisher", "Smart App Control",
                     "Get-FileHash", "sha256sum", "coal-gui", "SKILL.md", "/issues", "GPL-3.0-only",
                     # the owner asked for no DOI and named authors, so the page cites Newton and Murati and points
                     # to CITATION.cff
                     "Newton C, Murati K (2026)", "CITATION.cff"):
            self.assertIn(must, self.all)  # install, guide and cite pages
        # the owner asked for no DOI or Zenodo mention and for the large sidebar logo to go
        for gone in ("Zenodo", "placeholder", "DOI", 'class="logo"', "coal-logo.png"):
            for page, html in self.pages.items():
                self.assertNotIn(gone, html, page)

    def test_accessibility_basics(self):
        for page, html in self.pages.items():  # every page
            self.assertRegex(html, r'<html lang="en">', page)
            self.assertEqual(len(re.findall(r"<h1[\s>]", html)), 1, page)
            self.assertIn('<main id="main">', html, page)
            self.assertIn('href="#main"', html, page)
            for img in re.findall(r"<img\b[^>]*>", html):
                if 'src="assets/coal-mark-dark.svg"' in img:  # decorative, next to the text "COAL" in the same link
                    self.assertIn('alt=""', img)
                    continue
                self.assertRegex(img, r'\salt="[^"]{10,}"', img)
            self.assertRegex(html, r'<a class="home" href="\./"><img src="assets/coal-mark-dark\.svg"[^>]*><span>COAL</span></a>', page)
        css = _read("site", "style.css")
        self.assertIn("prefers-reduced-motion", css)
        # the page is dark only (no light variant)
        self.assertIn("color-scheme: dark;", css)
        self.assertNotIn("prefers-color-scheme", re.sub(r"/\*.*?\*/", "", css, flags=re.S))
        self.assertIn(":focus-visible", css)

    def test_same_header_and_footer_on_every_page(self):
        """One header and footer for the whole site; only the current page's nav link is marked."""
        def chrome(html):
            head = re.search(r'<header class="top">.*?</header>', html, re.S).group(0)
            foot = re.search(r'<footer class="foot".*?</footer>', html, re.S).group(0)
            return head.replace(' aria-current="page"', ""), foot
        first = chrome(self.index)
        for page, html in self.pages.items():
            self.assertEqual(chrome(html), first, f"{page}: header or footer differs from index.html")
            current = re.findall(r'<a href="([^"]+)" aria-current="page">', html)
            self.assertEqual(current, [] if page in ("index.html", "differences.html") else [page], page)
            if page not in ("index.html", "differences.html"):  # "Why use it" is reached from the home page
                self.assertIn(f'href="{page}"', first[0], f"{page} is not in the header")

    def test_font_is_self_hosted_with_its_licence(self):
        """The only font file is ours, small, and the OFL travels with it."""
        css = _read("site", "style.css")
        self.assertEqual(re.findall(r'url\("([^"]+)"\)', css), list(FONTS))
        for f in FONTS:  # Latin subsets of Literata (regular, semibold, italic)
            self.assertLess(os.path.getsize(os.path.join(SITE, f)), 40 * 1024, f)
        self.assertIn("SIL Open Font License", _read("site", "assets", "fonts", "OFL.txt"))
        self.assertIn("Literata", _read("site", "assets", "fonts", "OFL.txt"))  # the licence travels with the font

    def test_size_budget(self):
        total = 0
        for d, _, files in os.walk(SITE):
            for f in files:
                if not f.endswith((".png", ".jpg", ".jpeg", ".webp", ".gif", ".md")):
                    total += os.path.getsize(os.path.join(d, f))
        self.assertLess(total, 200 * 1024, f"site is {total} bytes without images (budget 200 KB)")


def _load_script(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "scripts", name + ".py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _script():
    spec = importlib.util.spec_from_file_location("site_for_release", os.path.join(ROOT, "scripts", "site_for_release.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class TestSiteForRelease(unittest.TestCase):
    def setUp(self):
        self.m = _script()

    def test_older_release_hides_newer_features_and_shows_its_version(self):
        out = self.m.rewrite(_read("site", "install.html"), "0.4.1")
        self.assertNotIn("Connect COAL to Claude", out)
        self.assertNotIn("data-since", out)
        self.assertIn("<li>Start <strong>COAL</strong>", out)  # the next step is untouched
        idx = self.m.rewrite(_read("site", "index.html"), "0.4.1")
        self.assertIn('<span id="version" data-version="0.4.1">0.4.1</span>', idx)
        self.assertIn("<span data-version-ref>0.4.1</span>", self.m.rewrite(_read("site", "cite.html"), "0.4.1"))

    def test_same_or_newer_release_keeps_features(self):
        for v in ("0.4.2", "0.5.0", "1.0"):
            out = self.m.rewrite(_read("site", "install.html"), v)
            self.assertIn("Connect COAL to Claude", out)
            self.assertNotIn("data-since", out)

    def test_idempotent(self):
        once = self.m.rewrite(_read("site", "install.html"), "0.4.1")
        self.assertEqual(self.m.rewrite(once, "0.4.1"), once)

    def test_bad_arguments(self):
        self.assertEqual(self.m.main(["latest", "x.html"]), 2)
        self.assertEqual(self.m.main(["0.4.1"]), 2)

    def test_pages_workflow_runs_it_on_the_deployed_copy(self):
        with open(os.path.join(ROOT, ".github", "workflows", "pages.yml"), encoding="utf-8") as fh:
            wf = fh.read()
        self.assertIn('python scripts/site_for_release.py "${tag#v}" _site/*.html', wf)
        self.assertLess(wf.index("site_for_release"), wf.index("upload-pages-artifact"))


class TestSiteWithoutRelease(unittest.TestCase):
    def test_downloads_replaced_when_no_release(self):
        # every page with release links (index.html and install.html)
        m = _load_script("site_without_release")
        with_links = []
        for name in sorted(os.listdir(os.path.join(ROOT, "site"))):
            if not name.endswith(".html"):
                continue
            with open(os.path.join(ROOT, "site", name), encoding="utf-8") as fh:
                page = fh.read()
            if not m.needs_rewrite(page):
                continue
            with_links.append(name)
            out = m.rewrite(page)
            self.assertNotIn("releases/latest/download/", out, name)
            self.assertNotIn("releases/download/", out, name)
            self.assertIn("First release coming soon", out, name)
            self.assertIn('id="no-release"', out, name)  # the target of every replaced link
            self.assertEqual(out.count("<main"), 1, name)
            self.assertFalse(m.needs_rewrite(out), name)
            self.assertEqual(m.rewrite(out), out)  # idempotent
        self.assertIn("index.html", with_links)
        self.assertIn("install.html", with_links)
        with open(os.path.join(ROOT, "site", "index.html"), encoding="utf-8") as fh:
            self.assertIn('pip install "coal-align @ git+https://github.com/plumat-coal/COAL" "coal-align-app[gui] @ '
                          'git+https://github.com/plumat-coal/COAL#subdirectory=apps/coal-app"', m.rewrite(fh.read()))

    def test_script_rewrites_every_page_given(self):
        import shutil
        import tempfile
        m = _load_script("site_without_release")
        with tempfile.TemporaryDirectory() as t:
            for name in ("index.html", "install.html", "cite.html"):
                shutil.copy(os.path.join(ROOT, "site", name), t)
            with open(os.path.join(t, "cite.html"), encoding="utf-8") as fh:
                cite = fh.read()
            paths = [os.path.join(t, n) for n in ("index.html", "install.html", "cite.html")]
            self.assertEqual(m.main(paths), 0)
            first = []
            for p in paths:
                with open(p, encoding="utf-8") as fh:
                    first.append(fh.read())
            self.assertNotIn("releases/latest/download/", first[0] + first[1])
            self.assertEqual(first[2], cite)  # no download links: left exactly as it was
            self.assertEqual(m.main(paths), 0)
            for p, before in zip(paths, first):
                with open(p, encoding="utf-8") as fh:
                    self.assertEqual(fh.read(), before)  # second run changes nothing

    def test_pages_workflow_checks_for_a_release(self):
        with open(os.path.join(ROOT, ".github", "workflows", "pages.yml"), encoding="utf-8") as fh:
            wf = fh.read()
        self.assertIn("gh release view", wf)
        self.assertIn("scripts/site_without_release.py _site/*.html", wf)  # every page
        self.assertIn("GH_TOKEN: ${{ github.token }}", wf)
        with open(os.path.join(ROOT, "site", "README.md"), encoding="utf-8") as fh:
            self.assertIn("before merging to `main`", fh.read())


class TestCoalrInstallLines(unittest.TestCase):
    """The five places that tell people how to install coalr name a release (REVIEW-COAL-3 H1): the release file
    coalr_<v>.tar.gz (from 0.6.1, the first release with one) or the tag v<v>; an unpinned install_github() of
    main only as the development version. At a release commit every pinned version is that release's and the
    "until <v> is published" notes are gone; in a development tree they name the release it becomes or an older one
    (docs/development.md "Releasing COAL")."""

    FILES = ("docs/r.md", "README.md", "packaging/r/coalr/README.md", "docs/tutorials/04-r-coalr.md",
             "packaging/orthologr/README.md")
    TARBALL = re.compile(r"https://github\.com/Plumat-COAL/COAL/releases/download/v([\d.]+)/coalr_([\d.]+)\.tar\.gz")
    TAG = re.compile(r'install_github\("Plumat-COAL/COAL@v([\d.]+)"')
    MAIN = re.compile(r'install_github\("Plumat-COAL/COAL"')

    @staticmethod
    def key(v):
        return tuple(int(x) for x in v.split("."))

    def test_pinned_to_a_release(self):
        rv = _load("release_version")
        version = rv.source_versions()[0]
        release = not rv.DEV_RE.match(version) and not rv.PRE_RE.search(version)
        target = rv.base_version(version)
        for rel in self.FILES:
            text = _read(*rel.split("/"))
            tarballs = self.TARBALL.findall(text)
            tags = self.TAG.findall(text)
            self.assertTrue(tarballs or tags, f"{rel}: no install line pinned to a release")
            for tag, name in tarballs:
                self.assertEqual(tag, name, f"{rel}: release v{tag} with coalr_{name}.tar.gz")
                self.assertGreaterEqual(self.key(tag), (0, 6, 1), f"{rel}: releases before 0.6.1 have no coalr file")
            for v in [t for t, _ in tarballs] + tags:
                if release:
                    self.assertEqual(v, version, f"{rel}: {version} is a release, but an install line names v{v}")
                else:
                    self.assertLessEqual(self.key(v), self.key(target), f"{rel}: v{v} is newer than this tree")
            for m in self.MAIN.finditer(text):
                around = text[max(0, m.start() - 200):text.find("\n", m.end())]
                self.assertIn("development version", around.lower(),
                              f"{rel}: an unpinned install_github() of main must say it is the development version")
            if release:
                self.assertNotIn(f"until {version} is published", text.lower(),
                                 f"{rel}: {version} is published with this commit: remove the note")

    def test_the_rule_itself(self):
        text = ('install.packages("https://github.com/Plumat-COAL/COAL/releases/download/v0.6.1/coalr_0.6.1.tar.gz")\n'
                'remotes::install_github("Plumat-COAL/COAL@v0.6.0", subdir = "packaging/r/coalr")\n')
        self.assertEqual(self.TARBALL.findall(text), [("0.6.1", "0.6.1")])
        self.assertEqual(self.TAG.findall(text), ["0.6.0"])
        self.assertEqual(self.MAIN.findall(text), [])
        self.assertEqual(len(self.MAIN.findall('install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")')), 1)


if __name__ == "__main__":
    unittest.main()
