"""The download website (site/), the release/pages workflows behind it, and the scripts that adapt the deployed copy to the latest release (or to none yet)."""
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
        # one version for both packages: the app's pyproject and the coal-align version it pins
        app, pin = rv.app_versions()
        self.assertEqual(app, init, "apps/coal-app/pyproject.toml has another version than the engine")
        self.assertEqual(pin, init, "apps/coal-app/pyproject.toml must depend on coal-align==<the same version>")
        self.assertEqual(rv.check("v" + init), (init, False))
        for bad in ("v9.9.9", init, "v" + init + "x", "", "release-" + init):
            with self.assertRaises(ValueError, msg=bad):
                rv.check(bad)
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "release_version.py"), "v9.9.9"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 1)
        self.assertIn("bump all four", r.stderr)
        self.assertIn("apps/coal-app/pyproject.toml", r.stderr)

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

    def test_ci_caches_pip_and_browsers(self):
        wf = _workflow("ci.yml")
        for name, job in wf["jobs"].items():
            for st in job["steps"]:
                if st.get("uses", "").startswith("actions/setup-python@"):
                    self.assertEqual(st["with"].get("cache"), "pip", name)
                    self.assertIn("pyproject.toml", st["with"].get("cache-dependency-path", ""), name)
        cache = [st for st in wf["jobs"]["explainer"]["steps"] if st.get("uses", "").startswith("actions/cache@")]
        self.assertEqual(len(cache), 1)
        self.assertEqual(cache[0]["with"]["path"], wf["jobs"]["explainer"]["env"]["PLAYWRIGHT_BROWSERS_PATH"])
        self.assertIn("playwright-1.48.0", cache[0]["with"]["key"])  # same pin as the install step

    def test_actions_pinned_to_commits(self):
        # ci.yml too (it uses actions/cache)
        for fname in ("release.yml", "ci.yml") + (("pages.yml",) if os.path.exists(os.path.join(WORKFLOWS, "pages.yml"))
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
            self.assertTrue(name in made or name in ("coal_align-${V}-py3-none-any.whl", "coal_align-${V}.tar.gz")
                            and "python -m build" in made, f"no build step makes {name}")
        sums = _steps(wf, "release")["SHA-256 of every file"]["run"]
        self.assertIn("> SHA256SUMS.txt", sums)


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
SITE_FILES = PAGES + ("assets/fonts/AtkinsonHyperlegibleNext-wght.woff2", "assets/fonts/OFL.txt") + ("index.html", "style.css", "main.js", "sw.js", "favicon.svg", "explainer.html", "img/broken-gene.webp",
              "assets/coal-mark.svg", "assets/coal-mark-dark.svg", "assets/coal-logo.png", "assets/coal-logo@2x.png", "assets/coal-logo-dark.png",
              "assets/coal-logo-dark@2x.png", "assets/favicon-32.png", "assets/apple-touch-icon.png",
              "assets/icon-512.png", "assets/coal-train.svg")
TEXT_EXT = (".html", ".css", ".js", ".svg", ".txt", ".json", ".webmanifest")
# Not requests: the SVG namespace, and the local app's address shown as text for people to open.
ALLOWED_NON_REPO_URLS = ("http://www.w3.org/2000/svg", "http://127.0.0.1:7860")


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

    def test_unsigned_and_smart_app_control_explained(self):
        for must in ("Windows protected your PC", "More info", "Run anyway", "Unknown publisher", "Smart App Control",
                     "Get-FileHash", "sha256sum", "coal-gui", "SKILL.md", "/issues", "MIT",
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
            self.assertEqual(current, [] if page == "index.html" else [page], page)
            if page != "index.html":
                self.assertIn(f'href="{page}"', first[0], f"{page} is not in the header")

    def test_font_is_self_hosted_with_its_licence(self):
        """The only font file is ours, small, and the OFL travels with it."""
        css = _read("site", "style.css")
        self.assertEqual(re.findall(r'url\("([^"]+)"\)', css), ["assets/fonts/AtkinsonHyperlegibleNext-wght.woff2"])
        size = os.path.getsize(os.path.join(SITE, "assets", "fonts", "AtkinsonHyperlegibleNext-wght.woff2"))
        self.assertLess(size, 60 * 1024)
        self.assertIn("SIL Open Font License", _read("site", "assets", "fonts", "OFL.txt"))
        self.assertIn('href="assets/fonts/OFL.txt"', self.index)

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


if __name__ == "__main__":
    unittest.main()
