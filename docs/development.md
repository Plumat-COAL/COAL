# Development

## Tests

| Where | What | Run |
|---|---|---|
| `tests/` | the engine (`coal`, `pal2nal`, `coal-cpp`): every stored PAL2NAL v14 reference output, options, large runs, real data; and the repository checks (`test_docs_*.py`: docs numbers, website, release workflows, skill; `test_repo_ci.py`: CI) | `python -m unittest discover -s tests` |
| `apps/coal-app/tests/` | the app: GUI form and layout, viewer, reports, MCP tools, settings, hosted Space, packaging, launchers (Gradio / mcp tests skip without the extra) | `python -m unittest discover -s apps/coal-app/tests` |
| `apps/coal-app/desktop/` | the Windows wrapper (`test_wrapper.py`) and an end-to-end smoke test of the local app (`test_smoke.py`) | `cd apps/coal-app/desktop && python -m unittest test_wrapper test_smoke` |
| `tests/fuzz_vs_pal2nal.py` | random comparison with the real PAL2NAL v14 (needs Perl; `python scripts/fetch_pal2nal.py`) | `PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python tests/fuzz_vs_pal2nal.py` |

What each comparison checks, the counts and the list of deliberate differences: [tests/README.md](../tests/README.md).

## Two packages

| Package | Folder | What | Depends on |
|---|---|---|---|
| `coal-align` | the repository root (`pyproject.toml`, `src/coal`, `src/coal_compat`) | the engine: `coal`, `pal2nal`; the example inputs | nothing (`[fast]`: numpy) |
| `coal-align-app` | `apps/coal-app` (`src/coal_app`, `desktop/`, `space/`) | web GUI, MCP server, viewer, Windows app, Hugging Face Space | `coal-align==` the same version; `[gui]`: Gradio, `[mcp]`: mcp |

A development install of both: `python -m pip install -e . -e "apps/coal-app[gui]"` (or `[mcp]`, in another
environment). Pipelines install the engine alone. The app runs the engine as `python -m coal` and imports none of its
internals. Build both: `python -m build .` and `python -m build apps/coal-app`.


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

1. Set the same version in `src/coal/__init__.py`, `pyproject.toml`, `apps/coal-app/pyproject.toml` (its
   `version` and its `coal-align==` dependency) and the `id="version"` span in `site/index.html` (a test fails if
   they differ; `scripts/release_version.py` refuses the tag), optionally add a `## <version>` section to
   `CHANGELOG.md` for the release notes, and commit.
2. Tag that commit: `git tag vX.Y.Z && git push origin vX.Y.Z`.
3. `.github/workflows/release.yml` checks the tag against the source version, builds and checks both wheels and
   sdists (engine and app) and the Windows app (`apps/coal-app/desktop/build.ps1`), writes `SHA256SUMS.txt` and
   publishes the GitHub Release, with
   version-free copies of the downloads. **Wait until it has published, then merge to `main`**, which redeploys the
   website (https://plumat-coal.github.io/COAL/, from `site/`, via `.github/workflows/pages.yml`). The site links to
   `releases/latest/download/<version-free name>`, so it offers the new release at once.
4. **PyPI** (`coal-align`, `coal-align-app`): the `pypi` job of the same workflow uploads the versioned wheels and
   sdists by [trusted publishing](https://docs.pypi.org/trusted-publishers/), after the GitHub Release exists. It
   runs only when the repository variable `PYPI_PUBLISH` is `true`. One-time setup: on PyPI, create both projects
   (or add a pending publisher) with owner `Plumat-COAL`, repository `COAL`, workflow `release.yml`, environment
   `pypi`; on GitHub, create the `pypi` environment (Settings > Environments; optionally with required reviewers)
   and the variable `PYPI_PUBLISH=true` (Settings > Secrets and variables > Actions > Variables). No API token is
   stored anywhere. To publish a release that was built before the variable was set, re-run its workflow from the
   Actions tab (`workflow_dispatch` with the tag).

The Windows app is not code-signed yet (the release notes say so, with the SmartScreen and Smart App Control steps).
Details and the one-time setup checklist: `site/README.md`.


## Requirements

- Python 3.10+ (the Windows app needs nothing else)
- A C++17 compiler and CMake 3.16+, only for building `coal-cpp`
- Perl 5, only for running the PAL2NAL reference
