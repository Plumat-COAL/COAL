# COAL website and releases

<!-- WEB-5: new file. See overnight/integration/CHANGES.md -->

`site/` is the download page at **https://plumat-coal.github.io/COAL/**: plain HTML, CSS and JavaScript, no build
step, nothing loaded from other sites (no CDNs, fonts, trackers or analytics). About 170 KB in all, 67 KB without
images. <!-- CHANGED(WEB-8): was "About 50 KB without the screenshot". See overnight/integration/CHANGES.md -->

| File | What it is |
|---|---|
| `index.html` | the page |
| `style.css` | styles, light and dark (follows the system setting) |
| `main.js` | picks the main download for the visitor's system, copy buttons, registers `sw.js` |
| `sw.js` | keeps a copy of the site's own files so it opens offline after one visit (never touches downloads) |
| `explainer.html` | a copy of `pitch/explainer.html` (the test fails if they differ: copy it again after editing the original) |
| `img/broken-gene.webp` | screenshot of the alignment view (`pitch/backup/1_broken_gene_frameshift.png`, lossless WebP) |
| `favicon.svg` | tab icon: one 3 x 3 codon block from the logo (lighter charcoal in a dark browser theme) |
| `assets/coal-logo*.png` | the owner's logo, 176 px wide at 1x and 2x; `-dark` has the wordmark recoloured to #ECE8E3 |
| `assets/coal-mark.svg` | the logo's squares as SVG (written by `desktop/make_icon.py --site site`) |
| `assets/coal-train.svg` | the owner's train drawing as SVG; the page inlines a copy in the footer so CSS can recolour it in dark mode |
| `assets/favicon-32.png`, `apple-touch-icon.png`, `icon-512.png` | PNG icons (written by `desktop/make_icon.py --site site`) |

<!-- CHANGED(WEB-8): logo, mark, train and icon files added; screenshot now WebP. See overnight/integration/CHANGES.md -->

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

1. Bump the version in `src/coal/__init__.py`, `pyproject.toml` and `index.html` (they must agree), optionally add a
   `## <version>` section to `CHANGELOG.md`, commit, and merge to `main`.
2. Tag and push: `git tag v0.3.0 && git push origin v0.3.0`.
3. The Release workflow builds and checks the wheel, sdist, Windows installer and portable zip, writes
   `SHA256SUMS.txt` and publishes the GitHub Release. The website's buttons point at it as soon as it is published.

A tag that doesn't match the source version fails at the first step. Tags like `v0.4.0rc1` become pre-releases,
which `releases/latest` (and so the website) ignores. To rebuild without publishing, run the workflow by hand
(Actions > Release > Run workflow) with the tag and "publish" off: the files are kept as workflow artifacts.

## The website deploy

`.github/workflows/pages.yml` runs the site tests and deploys `site/` (without this README) on every push to
`main` that touches the site, and can be run by hand.

## Preview locally

    python -m http.server -d site 8000      # then open http://127.0.0.1:8000/

To see the OS-specific variants, use "Show the download for Windows, macOS or Linux" under the button.

## One-time setup (repository admin)

- [ ] **Enable Pages:** Settings > Pages > Build and deployment > Source: **GitHub Actions**. Then run Actions >
      Pages > Run workflow once (or push to `main`).
- [ ] **First release:** push a tag `v<version>` matching `src/coal/__init__.py` (now `v0.3.0`). Until a release
      exists, every download link on the site gives a 404.
- [ ] Check the first release run's Windows job log (the build has never run on GitHub) and download the installer
      once to see the SmartScreen screens match the drawing on the page.
- [ ] Optional, later: **custom domain**: Settings > Pages > Custom domain, plus a DNS `CNAME` record. Nothing in
      the site needs changing (links are relative).
- [ ] Optional: **Zenodo DOI**: connect the repository at Zenodo (GitHub integration), make a release, then replace
      the placeholder in the "How to cite" section of `index.html` and add a `CITATION.cff`.
- [ ] Optional: a code-signing certificate would remove the SmartScreen warning; then update the Windows steps and
      the release notes text in `scripts/release_notes.py`.
