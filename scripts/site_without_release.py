"""Rewrite the assembled site for a deploy before the first release exists (FIN-A9).

    python scripts/site_without_release.py _site/index.html

Every `releases/latest/download/...` link would give a 404 until the first GitHub Release is published, so
.github/workflows/pages.yml runs this on the copy it deploys (never on site/ itself) when `gh release view` finds no
published release: download links become a "first release coming soon" note, the pip commands install from the
repository instead, and a short notice goes at the top of the page. Running it twice changes nothing more.
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
    # CHANGED(WEB-11): links first, keeping their text ("Download source code: first release coming soon"); was the
    # sdist URL first, which turned the source-code button into a broken href="git+https://..." link.
    # See overnight/integration/CHANGES.md
    html = LINK.sub(lambda m: f'<a class="pending" href="#no-release">{m.group(1)}: {NOTE.lower()}</a>', html)
    html = re.sub(re.escape(REPO) + r"/releases/(?:latest/)?download/coal_src\.tar\.gz", "git+" + REPO, html)
    if 'id="no-release"' not in html:
        html = re.sub(r'(<main\b[^>]*>)', lambda m: m.group(1) + "\n" + BANNER, html, count=1)
    return html


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[2].strip(), file=sys.stderr)
        return 2
    with open(argv[0], encoding="utf-8") as fh:
        html = fh.read()
    with open(argv[0], "w", encoding="utf-8", newline="\n") as fh:
        fh.write(rewrite(html))
    print(f"{argv[0]}: download links replaced (no published release yet)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
