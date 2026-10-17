# COAL website and releases

<!-- WEB-5: new file. See overnight/integration/CHANGES.md -->

`site/` is the website at **https://plumat-coal.github.io/COAL/**: plain HTML, CSS and JavaScript, no build
step, nothing loaded from other sites (no CDNs, font services, trackers or analytics). About 225 KB in all.
<!-- CHANGED(WEB-13): was one page, "About 170 KB in all, 67 KB without images" --> <!-- CHANGED(WEB-8): was "About 50 KB without the screenshot". See overnight/integration/CHANGES.md -->

| File | What it is |
|---|---|
| `index.html` | home: what COAL is, the three downloads, a short example, three reasons, a cite pointer |
| `install.html` | install steps for Windows (SmartScreen, Smart App Control), macOS, Linux, command line only; files and checksums |
| `guide.html` | example run, warnings and errors, limitations, links to the README and help |
| `differences.html` | differences from PAL2NAL v14 (the test counts) and the speed table |
| `cite.html` | how to cite COAL and PAL2NAL (the version is in a `data-version-ref` span) |
| `style.css` | styles for every page; dark only, light when printed |
| `assets/fonts/` | Atkinson Hyperlegible Next, variable weight 400 to 800, Latin subset as WOFF2 (25 KB), with its licence `OFL.txt` (SIL Open Font License 1.1) |
| `main.js` | copy buttons, registers `sw.js`; all three downloads are always shown (WEB-11), and CHANGED(WEB-13): the card for the visitor's system gets a "For this computer" mark, and old one-page links (`#windows`, `#speed`, ...) are sent to the page that now holds that section |
| `sw.js` | keeps a copy of the site's own files so it opens offline after one visit (never touches downloads) |
| `explainer.html` | a copy of `pitch/explainer.html` (the test fails if they differ: copy it again after editing the original) |
| `img/broken-gene.webp` | screenshot of the alignment view (`pitch/backup/1_broken_gene_frameshift.png`, lossless WebP) |
| `favicon.svg` | tab icon: one 3 x 3 codon block from the logo (lighter charcoal in a dark browser theme) |
| `assets/coal-logo*.png` | the owner's logo, 176 px wide at 1x and 2x; `-dark` has the wordmark recoloured to #ECE8E3. Not shown on the page since WEB-11 (kept for previews and other uses) |
| `assets/coal-mark.svg` | the logo's squares as SVG (written by `desktop/make_icon.py --site site`); the small mark in the page header |
| `assets/coal-train.svg` | the owner's train drawing as SVG; the page inlines a copy in the footer, beside the authors' names, so CSS can recolour it in dark mode |
| `assets/favicon-32.png`, `apple-touch-icon.png`, `icon-512.png` | PNG icons (written by `desktop/make_icon.py --site site`) |

<!-- CHANGED(WEB-8): logo, mark, train and icon files added; screenshot now WebP. See overnight/integration/CHANGES.md -->
<!-- CHANGED(WEB-13): several pages with one shared header and footer (a test checks they are identical; when you
     change the header or footer, change it in every page). The font was made from google/fonts
     ofl/atkinsonhyperlegiblenext with fontTools: `fontTools.varLib.instancer ... wght=400:800`, then `pyftsubset`
     to Basic Latin, Latin-1 and common punctuation, `--flavor=woff2`. -->

All internal links are relative, so the same files work at `/COAL/` and later at a custom domain.

## How downloads stay current

The download links never name a version. The page shows the current version once (`<span id="version">`), and
`tests_app/test_site.py` fails when it differs from `src/coal/__init__.py`, so bump it together with the package.
<!-- CHANGED(WEB-8): was "The page never names a version." See overnight/integration/CHANGES.md -->
Every download link has the form

    https://github.com/plumat-coal/COAL/releases/latest/download/<file>

which GitHub redirects to that file in the newest release (pre-releases excluded). The release workflow
(`.github/workflows/release.yml`) therefore uploads version-free copies next to the versioned files:

| Version-free (linked from the site) | Versioned (same content) |
|---|---|
| `COAL-Setup-windows-x64.exe` | `COAL-Setup-<version>.exe` |
| `COAL-windows-x64.zip` | `COAL-<version>-windows-x64.zip` |
| `coal_src.tar.gz` | `coal-<version>.tar.gz` |
| `SHA256SUMS.txt` | (covers every file) |
| (none) | `coal-<version>-py3-none-any.whl` |

Two pip details behind those names (tested with pip 24.0 through a local stand-in for GitHub's redirect):

- A wheel can't have a version-free name: pip refuses `coal-py3-none-any.whl` as "not a valid wheel filename". So
  the site's pip command installs the **sdist** alias; the release notes give the exact versioned wheel URL.
- The sdist alias has **no hyphen** on purpose. pip caches the wheel it builds from an archive URL whose name looks
  like `name-version`, so `coal-latest.tar.gz` kept installing the old version after a new release;
  `coal_src.tar.gz` doesn't get cached, so re-running the command upgrades.

`tests_app/test_site.py` checks that every `latest/download/` link on the page is in the release workflow's file
list, that the page has no other external URLs, that internal links resolve, and that every measured number on the
page (`data-check`) is quoted from `README.md`.

## Releasing

<!-- CHANGED(FIN-A9): order changed: tag and let the Release run publish BEFORE merging to main (merging deploys the
site at once). See overnight/integration/CHANGES.md -->
1. Bump the version in `src/coal/__init__.py`, `pyproject.toml` and `index.html` (they must agree; also the
   citation's `data-version-ref` span in `cite.html`, and `version:` in `CITATION.cff`), optionally add a
   `## <version>` section to `CHANGELOG.md`, and commit.
2. Tag that commit and push the tag: `git tag v0.4.0 && git push origin v0.4.0`.
3. The Release workflow builds and checks the wheel, sdist, Windows installer and portable zip, writes
   `SHA256SUMS.txt` and publishes the GitHub Release. **Wait until it has published, before merging to `main`**:
   merging deploys the website at once, and its buttons point at the newest published release.
4. Merge to `main`. If no published release exists when the site is deployed, `pages.yml` replaces the download
   buttons with a "First release coming soon" note and a pip install from the repository
   (`scripts/site_without_release.py`), so visitors never get a 404.

A tag that doesn't match the source version fails at the first step. Tags like `v0.4.0rc1` become pre-releases,
which `releases/latest` (and so the website) ignores. To rebuild without publishing, run the workflow by hand
(Actions > Release > Run workflow) with the tag and "publish" off: the files are kept as workflow artifacts.

## The website deploy

`.github/workflows/pages.yml` runs the site tests and deploys `site/` (without this README) on every push to
`main` that touches the site, and can be run by hand.

## Preview locally

    python -m http.server -d site 8000      # then open http://127.0.0.1:8000/

To see the site as deployed before the first release: copy `site/` somewhere, run
`python scripts/site_without_release.py <copy>/*.html` and serve the copy. <!-- CHANGED(WEB-13): was <copy>/index.html -->

## One-time setup (repository admin)

- [ ] **Enable Pages:** Settings > Pages > Build and deployment > Source: **GitHub Actions**. Then run Actions >
      Pages > Run workflow once (or push to `main`).
- [ ] **First release, before merging to `main`:** push a tag `v<version>` matching `src/coal/__init__.py` (now
      `v0.4.0`) and wait for the Release run to publish. Until a release exists, the deployed site shows "First
      release coming soon" instead of download buttons (FIN-A9).
- [ ] Check the first release run's Windows job log (the build has never run on GitHub) and download the installer
      once to see the SmartScreen screens match the drawing on the page.
- [ ] Optional, later: **custom domain**: Settings > Pages > Custom domain, plus a DNS `CNAME` record. Nothing in
      the site needs changing (links are relative).
<!-- CHANGED(WEB-11): the Zenodo DOI item is gone: the owner decided on no DOI; the page cites the authors, version
     and repository, and CITATION.cff is at the repository root -->
- [ ] Optional: a code-signing certificate would remove the SmartScreen warning; then update the Windows steps and
      the release notes text in `scripts/release_notes.py`.
