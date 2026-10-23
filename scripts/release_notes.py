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
# the engine wheels with coal-cpp inside that release.yml builds (job cpp-wheels), in its order
PLATFORM_TAGS = ("manylinux_2_28_x86_64", "manylinux_2_28_aarch64", "macosx_11_0_universal2", "win_amd64")
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
    app_wheel = f"coal_align_app-{version}-py3-none-any.whl"
    platform_wheels = ", ".join(f"`coal_align-{version}-py3-none-{t}.whl`" for t in PLATFORM_TAGS)
    return f"""COAL {version}: turns a protein alignment plus the matching coding DNA into a codon alignment
(for dN/dS and DNA phylogenies). A fast, MIT-licensed replacement for PAL2NAL.

Website and install steps: https://plumat-coal.github.io/COAL/

## Downloads

| File | What it is |
|---|---|
| `COAL-Setup-{version}.exe` (same file: `COAL-Setup-windows-x64.exe`) | Windows 10/11 installer, 64-bit, installs for you only (no administrator rights) |
| `COAL-{version}-windows-x64.zip` (same file: `COAL-windows-x64.zip`) | Windows portable app: unzip anywhere, run `COAL\\COAL.exe` |
| {platform_wheels} | Python wheels of the engine, `coal-align`, with the C++ engine `coal-cpp` inside (Linux x86_64 and aarch64 with glibc 2.28+, macOS 11+ on Apple silicon and Intel, Windows x64; Python 3.10+; no dependencies). `pip install` picks the one for your platform |
| `{wheel}` | Python wheel of the engine without `coal-cpp` (any platform with Python 3.10+; the fallback pip uses where no wheel above fits) |
| `coal_align-{version}.tar.gz` (same file: `coal_src.tar.gz`) | Python source distribution of the engine |
| `{app_wheel}` | Python wheel of the app, `coal-align-app` (web GUI and MCP server; needs the engine wheel) |
| `coal_align_app-{version}.tar.gz` (same file: `coal_app_src.tar.gz`) | Python source distribution of the app |
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

Needs Python 3.10 or newer. In a virtual environment, the engine and the app:

```
python -m pip install "coal-align @ {REPO}/releases/download/v{version}/coal_align-{version}.tar.gz" \\
    "coal-align-app[gui] @ {REPO}/releases/download/v{version}/{app_wheel}"
coal-gui        # the app, at http://127.0.0.1:7860
coal --help     # the command line
```

For the command line alone (pipelines), install only the first package, `coal-align`: it has no dependencies.
To get the C++ engine with it, install the wheel for your platform from the table above instead of the
`.tar.gz` (download it, then `python -m pip install coal_align-{version}-py3-none-<platform>.whl`); `coal --cpp-path`
then prints where `coal-cpp` is. On PyPI, `python -m pip install coal-align` chooses that wheel by itself.
A plain `pip install coal` installs a different project.

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
