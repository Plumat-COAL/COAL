"""Write the GitHub Release notes for one version (used by .github/workflows/release.yml).

    python scripts/release_notes.py 0.4.0 notes.md

Notes = a fixed part (what to download, "not code-signed" with the SmartScreen / Smart App Control caveat, the pip
alternative, how to check SHA256SUMS.txt) + the "## <version>" section of CHANGELOG.md if there is one.
Exit status 0 when a CHANGELOG section was found, 3 when not (the workflow then also asks GitHub to generate
notes from the merged pull requests). Any other status is an error.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = "https://github.com/plumat-coal/COAL"
NO_CHANGELOG = 3


def changelog_section(version, path=os.path.join(ROOT, "CHANGELOG.md")):
    """The body under a heading like '## 0.3.0', '## v0.3.0', '## [0.3.0] - 2026-10-18'; None if absent."""
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    head = re.compile(r"(?m)^##\s+\[?v?" + re.escape(version) + r"\]?(?:\s.*)?$")
    m = head.search(text)
    if not m:
        return None
    nxt = re.compile(r"(?m)^##\s").search(text, m.end())
    body = text[m.end():nxt.start() if nxt else len(text)].strip()
    return body or None


def fixed_part(version):
    wheel = f"coal_align-{version}-py3-none-any.whl"
    return f"""COAL {version}: turns a protein alignment plus the matching coding DNA into a codon alignment
(for dN/dS and DNA phylogenies). A fast, MIT-licensed replacement for PAL2NAL.

Website and install steps: https://plumat-coal.github.io/COAL/

## Downloads

| File | What it is |
|---|---|
| `COAL-Setup-{version}.exe` (same file: `COAL-Setup-windows-x64.exe`) | Windows 10/11 installer, 64-bit, installs for you only (no administrator rights) |
| `COAL-{version}-windows-x64.zip` (same file: `COAL-windows-x64.zip`) | Windows portable app: unzip anywhere, run `COAL\\COAL.exe` |
| `{wheel}` | Python wheel (Windows, macOS, Linux; Python 3.10+) |
| `coal_align-{version}.tar.gz` (same file: `coal_src.tar.gz`) | Python source distribution |
| `SHA256SUMS.txt` | SHA-256 checksums of every file above |

The version-free names always point at the newest release through
`{REPO}/releases/latest/download/<name>`, which is what the website links to.

## The Windows app is not code-signed

COAL has no code-signing certificate yet, so Windows doesn't know the publisher:

- **Microsoft Defender SmartScreen** shows "Windows protected your PC". Click **More info**, check that it says
  `COAL-Setup-...exe` and "Publisher: Unknown publisher", then **Run anyway**. Edge may first say the file
  "isn't commonly downloaded": choose **...** > **Keep** > **Show more** > **Keep anyway**.
- **Smart App Control** (Windows 11) blocks unsigned apps outright, with no "Run anyway" and no per-app exception.
  If it is on, use pip instead (below); don't turn Smart App Control off for COAL.

## Install with pip instead (Windows, macOS, Linux)

Needs Python 3.10 or newer. In a virtual environment:

```
python -m pip install "coal-align[gui] @ {REPO}/releases/download/v{version}/{wheel}"
coal-gui        # the app, at http://127.0.0.1:7860
coal --help     # the command line
```

COAL's package name is `coal-align`: a plain `pip install coal` installs a different project.

## Check a download

Compare with the line for that file in `SHA256SUMS.txt`:

- Windows (PowerShell): `Get-FileHash .\\COAL-Setup-windows-x64.exe -Algorithm SHA256`
- Linux: `sha256sum -c --ignore-missing SHA256SUMS.txt`
- macOS: `shasum -a 256 coal_src.tar.gz`
"""


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    version, out = argv
    section = changelog_section(version)
    text = fixed_part(version)
    if section:
        text += f"\n## Changes in {version}\n\n{section}\n"
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    return 0 if section else NO_CHANGELOG


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
