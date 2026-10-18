"""Make the deployed site describe the release people actually download (WEB-18).

    python scripts/site_for_release.py 0.4.1 _site/*.html

The site is built from main, but its download buttons give the newest published release. Between a merge to main and
the release (or if a release build fails), the two differ: the site said "Version 0.4.2" and described the 0.4.2
installer's "Connect COAL to Claude Desktop" box while the button gave 0.4.1. .github/workflows/pages.yml runs this on
the copy it deploys (never on site/ itself) with the latest release's version:

- the version shown (`<span id="version" data-version=...>`, `<span data-version-ref>`) becomes that version;
- an element marked `data-since="X.Y.Z"` (an `<li>`, `<p>`, `<div>` or `<section>`, not nested in one of the same
  kind) is removed when the release is older than X.Y.Z, so the page never describes a feature the download lacks.

Running it twice changes nothing more. With the same version as main, only the data-since markers go.
Rollback: delete the "Match the site to the latest release" step in pages.yml (or revert the WEB-18 commit).
"""
from __future__ import annotations

import re
import sys

VERSION = re.compile(r"^\d+(?:\.\d+)*$")


def key(v: str) -> tuple[int, ...]:
    return tuple(int(x) for x in v.split("."))


def rewrite(html: str, version: str) -> str:
    html = re.sub(r'(<span id="version" data-version=")[^"]*(">)[^<]*(</span>)',
                  lambda m: f"{m.group(1)}{version}{m.group(2)}{version}{m.group(3)}", html)
    html = re.sub(r"(<span data-version-ref>)[^<]*(</span>)", lambda m: f"{m.group(1)}{version}{m.group(2)}", html)

    def since(m: re.Match) -> str:
        return "" if key(version) < key(m.group(2)) else m.group(0).replace(f' data-since="{m.group(2)}"', "", 1)

    return re.sub(r'[ \t]*<(li|p|div|section)\b[^>]*?\sdata-since="([\d.]+)"[^>]*>.*?</\1>(?:[ \t]*\n)?', since, html,
                  flags=re.S)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2 or not VERSION.match(argv[0].lstrip("v")):
        print("usage: python scripts/site_for_release.py VERSION PAGE.html ...", file=sys.stderr)
        return 2
    version = argv[0].lstrip("v")
    for path in argv[1:]:
        with open(path, encoding="utf-8") as fh:
            html = fh.read()
        new = rewrite(html, version)
        if new != html:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(new)
            print(f"{path}: matched to release {version}")
        else:
            print(f"{path}: unchanged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
