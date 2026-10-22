"""Rewrite the assembled site for a deploy before the first release exists (FIN-A9).

    python scripts/site_without_release.py _site/*.html

Every `releases/latest/download/...` link would give a 404 until the first GitHub Release is published, so
.github/workflows/pages.yml runs this on the copy it deploys (never on site/ itself) when `gh release view` finds no
published release: download links become a "first release coming soon" note, the pip commands install from the
repository instead, and a short notice goes at the top of the page. Running it twice changes nothing more.
Takes any number of pages; a page with no release download link
(the guide, the cite page, the explainer) is left exactly as it is.
"""
from __future__ import annotations

import re
import sys

REPO = "https://github.com/plumat-coal/COAL"
PIP = f'pip install "coal[gui] @ git+{REPO}"'
NOTE = "First release coming soon"
BANNER = (f'<p class="note" role="status" id="no-release">{NOTE}: there is nothing to download yet. Install from '
          f'source instead (Python 3.10 or newer): <code>{PIP}</code></p>')
LINK = re.compile(r'<a\b[^>]*href="' + re.escape(REPO) + r'/releases/(?:latest/)?download/[^"]*"[^>]*>(.*?)</a>', re.S)


def rewrite(html: str) -> str:
    # Links first, keeping their text ("Download source code: first release coming soon"); rewriting the sdist URL
    # first would turn the source-code button into a broken href="git+https://..." link.
    html = LINK.sub(lambda m: f'<a class="pending" href="#no-release">{m.group(1)}: {NOTE.lower()}</a>', html)
    html = re.sub(re.escape(REPO) + r"/releases/(?:latest/)?download/coal_src\.tar\.gz", "git+" + REPO, html)
    if 'id="no-release"' not in html:
        html = re.sub(r'(<main\b[^>]*>)', lambda m: m.group(1) + "\n" + BANNER, html, count=1)
    return html


def needs_rewrite(html: str) -> bool:
    """WEB-13: only pages that link to a release file get the notice."""
    return LINK.search(html) is not None or REPO + "/releases/latest/download/" in html


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__.strip().splitlines()[2].strip(), file=sys.stderr)
        return 2
    for path in argv:
        with open(path, encoding="utf-8") as fh:
            html = fh.read()
        if not needs_rewrite(html):
            print(f"{path}: no release download links, unchanged")
            continue
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(rewrite(html))
        print(f"{path}: download links replaced (no published release yet)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
