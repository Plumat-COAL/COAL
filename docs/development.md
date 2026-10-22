# Development

## Tests

| Where | What | Run |
|---|---|---|
| `tests/` | the engine (`coal`, `pal2nal`, `coal-cpp`): every stored PAL2NAL v14 reference output, options, large runs, real data | `python -m unittest discover -s tests` |
| `tests_app/` | the app: GUI logic, MCP tools, reports, docs, website, packaging (Gradio / mcp tests skip without the extra) | `python -m unittest discover -s tests_app` |
| `desktop/` | the Windows wrapper (`test_wrapper.py`) and an end-to-end smoke test of the local app (`test_smoke.py`) | `cd desktop && python -m unittest test_wrapper test_smoke` |
| `tests/fuzz_vs_pal2nal.py` | random comparison with the real PAL2NAL v14 (needs Perl; `python scripts/fetch_pal2nal.py`) | `PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python tests/fuzz_vs_pal2nal.py` |

What each comparison checks, the counts and the list of deliberate differences: [tests/README.md](../tests/README.md).


## Reference implementation (PAL2NAL)

COAL is an independent implementation. PAL2NAL v14 is used only as a black-box reference for differential testing
and benchmarking. PAL2NAL is distributed under the GNU General Public License v2.0 (GPL-2.0). COAL contains no
PAL2NAL code, so **PAL2NAL is not included or redistributed in this repository**, and its code and test data must
not be copied here.

To fetch it locally (into `PAL2NAL/`, which is git-ignored):

```
python scripts/fetch_pal2nal.py
```

If the official site can't be reached, the script falls back to the copy of `pal2nal.pl` in the public orthologr
R package on GitHub, accepted only if its SHA-256 equals that of the official script. If both downloads are blocked,
put `pal2nal.v14.tar.gz` (downloaded from the official site in a browser) in the repo root and run
`python scripts/fetch_pal2nal.py --archive pal2nal.v14.tar.gz`. Run it (requires Perl 5):

```
perl PAL2NAL/pal2nal.v14/pal2nal.pl pep.aln nuc.fasta -output fasta > out.fa 2> out.err
```

Reference: Suyama M, Torrents D, Bork P (2006) PAL2NAL: robust conversion of protein sequence alignments into the
corresponding codon alignments. *Nucleic Acids Res.* 34:W609–W612.

## Releasing COAL

1. Set the same version in `src/coal/__init__.py`, `pyproject.toml` and the `id="version"` span in
   `site/index.html` (a test fails if they differ), optionally add a `## <version>` section to `CHANGELOG.md` for
   the release notes, and commit.
2. Tag that commit: `git tag vX.Y.Z && git push origin vX.Y.Z`.
3. `.github/workflows/release.yml` checks the tag against the source version, builds and checks the wheel, the sdist
   and the Windows app (`desktop/build.ps1`), writes `SHA256SUMS.txt` and publishes the GitHub Release, with
   version-free copies of the downloads. **Wait until it has published, then merge to `main`**, which redeploys the
   website (https://plumat-coal.github.io/COAL/, from `site/`, via `.github/workflows/pages.yml`). The site links to
   `releases/latest/download/<version-free name>`, so it offers the new release at once.

The Windows app is not code-signed yet (the release notes say so, with the SmartScreen and Smart App Control steps).
Details and the one-time setup checklist: `site/README.md`.


## Requirements

- Python 3.10+ (the Windows app needs nothing else)
- A C++17 compiler and CMake 3.16+, only for building `coal-cpp`
- Perl 5, only for running the PAL2NAL reference
