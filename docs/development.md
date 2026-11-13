# Development

## Tests

| Where | What | Run |
|---|---|---|
| `tests/` | the engine (`coal`, `pal2nal`, `coal-cpp`): every stored PAL2NAL v14 reference output, options, large runs, real data; and the repository checks (`test_docs_*.py`: docs numbers, website, release workflows, skill; `test_repo_ci.py`: CI) | `python -m unittest discover -s tests` |
| `apps/coal-app/tests/` | the app: GUI form and layout, viewer, reports, MCP tools, settings, hosted Space, packaging, launchers (Gradio / mcp tests skip without the extra) | `python -m unittest discover -s apps/coal-app/tests` |
| `apps/coal-app/desktop/` | the Windows wrapper (`test_wrapper.py`) and an end-to-end smoke test of the local app (`test_smoke.py`) | `cd apps/coal-app/desktop && python -m unittest test_wrapper test_smoke` |
| `tests/fuzz_vs_pal2nal.py` | random comparison with the real PAL2NAL v14 (needs Perl; `python scripts/fetch_pal2nal.py`) | `PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python tests/fuzz_vs_pal2nal.py` |
| `tests/tools/engine_fuzz.py` | random `coal-cpp` vs Python engine comparison on malformed inputs and every option (bytes, exit code, report); run nightly with the PAL2NAL fuzz for both engines and an ASan/UBSan build by `.github/workflows/nightly.yml` (see [CONTRIBUTING.md](../CONTRIBUTING.md#nightly-checks)) | `COAL_CPP_BIN=build/coal-cpp python tests/tools/engine_fuzz.py --n 500 --seed 1` |

What each comparison checks, the counts and the list of deliberate differences: [tests/README.md](../tests/README.md).
The short practical guide (setup, the rule that both engines stay byte-identical, how to add a deliberate
difference, help-text regeneration): [CONTRIBUTING.md](../CONTRIBUTING.md). The public Python API and the 1.0
versioning promise: [python-api.md](python-api.md); the run report's fields and JSON Schema:
[report-schema.md](report-schema.md) (`tests/test_report_schema.py` validates both engines' reports against it).

## Test hooks

Some tests have to stop a run at a chosen point (to kill it there), or to make several runs reach the output lock
in the same instant. They do that with environment variables named `COAL_TEST_*`. None of them is part of COAL's
interface, and none does anything in an ordinary run.

**`coal-cpp`.** The hooks are compiled in only when the build is configured with `-DCOAL_TEST_HOOKS=ON` (a CMake
option, `OFF` by default). Release builds, the wheels' binaries, the Bioconda build and the desktop app's binary
are built without: their `coal-cpp` does not contain the code and never reads the variables. The test builds of CI
(`ci.yml` job `cpp`, the sanitizer job of `nightly.yml`) and a local build for the test suite switch the option on:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOAL_TEST_HOOKS=ON && cmake --build build --config Release
```

Even in such a build a hook acts only while `COAL_TEST_HOOKS=1` is set in the environment as well.

| Variable | Effect (test builds of `coal-cpp`) | Python engine |
| --- | --- | --- |
| `COAL_TEST_HOOKS=1` | switches the hooks below on | the same variable (`pipeline._test_hooks`) |
| `COAL_TEST_HOOKS=probe` | the binary prints `coal-cpp: test hooks compiled in` and exits 0: how the tests ask a binary | none |
| `COAL_TEST_PAUSE_AFTER_CHECKPOINT=k` | sleeps after the k-th checkpoint | the same variable |
| `COAL_TEST_PAUSE_BEFORE_CHECKPOINT=1` | sleeps before the first checkpoint, partial files written | none |
| `COAL_TEST_PAUSE_IN_REPORT=1` | sleeps while the report's temporary file exists | the same variable (`report.py`) |
| `COAL_TEST_LOCK_START_AT=<ms since the epoch>` | waits until then before it takes the output lock | none (the tests pass the time to their own starter) |
| `COAL_TEST_HOLD_LOCK_UNTIL=<path>` | holds the output lock until that file exists, then checks the lock is still its own | none |
| `COAL_TEST_LOCK_BEFORE_TAKEOVER_UNTIL=<file>` | between judging a lock stale and replacing it: creates `<file>.waiting`, then waits until `<file>` exists (at most 60 s), so a test can act inside the window without a fixed sleep | |
| `COAL_TEST_LOCK_BEFORE_TAKEOVER_MS=<ms>` | waits between judging a lock stale and replacing it | `pipeline._before_takeover`, `pipeline._in_takeover` (attributes a test sets) |
| `COAL_TEST_UNLINK_FAILS=N` | the first N removals of lock, guard and temporary files fail as on Windows while another process has the file open | `pipeline._unlink_fails` (an attribute a test sets) |
| `COAL_TEST_DNDS_FAIL_AFTER_GROUPS=k` | `coal-cpp codon-dnds` fails with a write error when k groups have been written (what it wrote so far is flushed first): a stage that dies in the middle, for `tests/test_pipeline_atomic.py` | none |

**Python engine.** Python has no build step to leave the hooks out of, so the three it reads from the environment
(`COAL_TEST_PAUSE_AFTER_CHECKPOINT` in `pipeline.py`, `COAL_TEST_PAUSE_IN_REPORT` in `report.py`, both only with
`COAL_TEST_HOOKS=1`) stay in the installed package: test-only, not an interface, and they may change without
notice. The other Python hooks are module attributes that only a test in the same process can set.

**Variables that are not hooks.** `COAL_TEST_DIAMOND` (the DIAMOND binary the tests use), `COAL_TEST_SLOW_FACTOR`
(longer time limits in the app's tests) and `COAL_LOCK_RACE_TRIALS` are read by the tests only, not by COAL.

**Tests and binaries without hooks.** A test that needs a hook asks the binary first (`tests/cpp_hooks.py`). With a
`coal-cpp` that has none it is skipped with that reason; when `COAL_REQUIRE_CPP` is set (CI) it fails instead, so
a test build without the option cannot pass unnoticed. `COAL_ALLOW_NO_TEST_HOOKS=1` turns that failure back into
a skip, for running the suite against a release binary on purpose.

## Two packages

| Package | Folder | What | Depends on |
|---|---|---|---|
| `coal-align` | the repository root (`pyproject.toml`, `src/coal`, `src/coal_compat`) | the engine: `coal`, `pal2nal`; the example inputs | nothing (`[fast]`: numpy) |
| `coal-align-app` | `apps/coal-app` (`src/coal_app`, `desktop/`, `space/`) | web GUI, MCP server, viewer, Windows app, Hugging Face Space | `coal-align==` the same version; `[gui]`: Gradio, `[mcp]`: mcp |

A development install of both: `python -m pip install -e . -e "apps/coal-app[gui]"` (or `[mcp]`, in another
environment). Pipelines install the engine alone. The app runs the engine as `python -m coal` and imports none of its
internals. Build both: `python -m build .` and `python -m build apps/coal-app`.


## Reference implementation (PAL2NAL)

COAL's codon-alignment engine is an independent implementation (COAL itself is GPL-3.0-only: `LICENSE`, `NOTICE`).
PAL2NAL v14 is used only as a black-box reference for differential testing
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
