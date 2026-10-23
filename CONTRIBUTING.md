# Contributing to COAL

COAL is two engines that must behave identically (Python in `src/coal`, C++ in `cpp/`), a `pal2nal` drop-in, and an
app in `apps/coal-app`. This page is the short version of how to work on it; the details are in
[docs/development.md](docs/development.md) and [tests/README.md](tests/README.md).

## Set up

```
git clone https://github.com/Plumat-COAL/COAL && cd COAL
python -m venv .venv && . .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e . -e "apps/coal-app[gui]"      # the engine (coal-align) and the app (coal-align-app)
python -m pip install psutil                            # tests only (the worker-process tests)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release -j 4   # coal-cpp
export COAL_CPP_BIN=$PWD/build/coal-cpp                 # Windows: build\Release\coal-cpp.exe
```

Python 3.10+; a C++17 compiler and CMake 3.16+ for `coal-cpp`; Perl 5 only for the PAL2NAL comparison. The `[mcp]`
extra may not install next to `[gui]` (pins): use a second environment for it. `numpy` (optional, `[fast]`) and
`pyyaml` (the workflow checks in `tests/test_docs_*`) are worth installing in a development environment; `psutil`
is needed only by the tests that watch COAL's worker processes (`tests/test_robustness.py`,
`tests/test_engine_extras.py`), which skip without it. It is not a dependency of the package.

`COAL_CPP_BIN` tells the tests where the binary is (default `./build/coal-cpp`); without it the C++ tests skip.
`COAL_IMPL=cpp` makes the reference comparison run the C++ engine; at run time `coal` hands large inputs to
`coal-cpp` when it finds one (`COAL_CPP=<path>`, `COAL_ENGINE=python|cpp` to force).

## Run the tests

The same commands as `.github/workflows/ci.yml`:

| What | Command |
|---|---|
| Engine (Python): goldens, options, report, API, docs and repository checks | `python -m unittest discover -s tests` |
| Engine, C++ engine included (goldens with `COAL_IMPL=cpp`, C++ vs Python, CLI, readers, tables, thread pool) | `COAL_IMPL=cpp COAL_CPP_BIN=$PWD/build/coal-cpp python -m unittest -v tests.test_differential tests.test_cpp_vs_python tests.test_cpp_cli tests.test_cpp_readers tests.test_cpp_tables tests.test_cpp_pool` |
| Everything in `tests/` with the binary built (what to run before a PR) | `COAL_CPP_BIN=$PWD/build/coal-cpp python -m unittest discover -s tests` |
| App | `python -m unittest discover -s apps/coal-app/tests` (`COAL_REQUIRE_EXTRAS=gui` or `mcp` makes the extra's tests fail instead of skip) |
| Windows wrapper and local-app smoke test | `cd apps/coal-app/desktop && python -m unittest test_wrapper test_smoke` |
| Random comparison with the real PAL2NAL v14 | `python scripts/fetch_pal2nal.py` then `PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl FUZZ_N=3000 python tests/fuzz_vs_pal2nal.py -v` |

`scripts/fetch_pal2nal.py` downloads PAL2NAL into `PAL2NAL/` (git-ignored; sha256-checked; never committed: PAL2NAL
is GPL-2.0 and COAL contains none of its code). CI runs the fuzz on every push (job `pal2nal-fuzz`, 3,000 cases) and
fails if it was skipped. `python .github/ci_unittest.py <unittest args>` is `python -m unittest` with GitHub
annotations; locally plain `unittest` is fine.

## The two engines stay byte-identical

`coal` and `coal-cpp` must write the same output bytes, the same warnings, the same run report (apart from
`runtime_seconds`) and the same exit codes for every input and option. This is enforced by:

- `tests/test_differential.py` with `COAL_IMPL=cpp`: all 285 stored PAL2NAL v14 reference checks through the C++
  engine;
- `tests/test_cpp_vs_python.py`: random inputs (all formats and options, report JSON included) through both engines,
  compared byte for byte (`COAL_RANDOM_N` changes the count);
- `tests/test_cpp_vs_reference.py`: the Python engine's own test modules (`test_decisions`, `test_report`, ...) run
  against a shim whose `coal` command is `coal-cpp`;
- `tests/test_cpp_cli.py`: the help text and the CLI's argument errors; `tests/test_report_schema.py`: both engines'
  reports against `docs/report.schema.json`.

So a change in behaviour is made **in both** `src/coal` and `cpp/src` in the same pull request (the C++ files carry a
"mirrors ..." comment naming the Python function they follow). If a change touches the command line (an option, its
help), regenerate the C++ help text from the Python parser and commit it:

```
python tests/test_cpp_cli.py --write-help          # rewrites cpp/src/help_text.hpp
```

## Changing behaviour vs PAL2NAL

COAL's output equals PAL2NAL v14's except for a documented list of deliberate differences, and the tests count them.
A new deliberate difference needs, in one pull request:

1. A reason other than "PAL2NAL does it that way": the current NCBI tables, a format's definition, PAL2NAL's own
   inconsistency, or a stated policy (see the head of [docs/differences.md](docs/differences.md)).
2. The behaviour in both engines, with a fixture if a stored golden can show it (`tests/fixtures/<case>/`, listed in
   `tests/cases.json`; `python scripts/make_golden.py` records PAL2NAL's output, which needs Perl and the fetched
   PAL2NAL; goldens are never edited by hand).
3. Its entry in the `DELIBERATE` table of `tests/test_differential.py` (`error`, `same_output`, `differs` or
   `nogap_code`, per case and variant) and, if the fuzz generator can produce it, a feature in the `DELIBERATE`
   dict of `tests/fuzz_vs_pal2nal.py` with a neutralised reference (the fuzz accepts a difference only when it is
   exactly the documented one).
4. The documentation: a bullet in `docs/differences.md`, a row in the "Deliberate differences" table of
   `tests/README.md`, and the updated counts in `tests/README.md` ("Test numbers": the printed summary line of
   `test_differential.py` is quoted there and checked by `tests/test_docs_numbers.py`, which also checks `README.md`,
   the website and the skill), plus a `CHANGELOG.md` line.

An *unexplained* difference found by the fuzz is a bug in COAL until shown otherwise.

## Pull requests

- One topic per PR; tests for the change (the examples in `docs/python-api.md` and `docs/report-schema.md` are
  tests too: `tests/test_python_api.py`, `tests/test_report_schema.py`).
- Run `COAL_CPP_BIN=$PWD/build/coal-cpp python -m unittest discover -s tests` and the app tests before pushing. CI
  runs them on Linux, Windows and macOS with Python 3.10 and 3.13, and on Linux with 3.14; the C++ build uses warnings as errors
  (`-DCOAL_WARNINGS_AS_ERRORS=ON`).
- Output bytes, exit codes, the report schema (`coal.report/1`: additive only) and the public Python API are
  covered by the versioning promise in [docs/python-api.md](docs/python-api.md#versioning-promise): a change to any
  of them is called out in the PR and in `CHANGELOG.md`.
- Documentation lives next to what it describes: `README.md` (the short tour and the Documentation table),
  `docs/usage.md` (every option), `docs/differences.md`, `docs/cpp-engine.md`, `tests/README.md` (what the tests
  check and the counts). Numbers quoted in docs are checked by `tests/test_docs_numbers.py` and `tests/test_docs_site.py`.
- Issue forms and the PR checklist are in `.github/`.

## Releasing

The version lives in four places and must be the same in all: `src/coal/__init__.py`, `pyproject.toml`,
`apps/coal-app/pyproject.toml` (its `version` and its `coal-align==` pin), plus the `id="version"` span in
`site/index.html`. `scripts/release_version.py vX.Y.Z` checks a tag against them (it is what
`.github/workflows/release.yml` runs; it refuses a tag that disagrees with the source) and prints the version:

```
python scripts/release_version.py v0.4.2       # prints 0.4.2, exit 0; exit 1 with the mismatch spelled out
```

Steps: bump the four places, add a `## <version>` section to `CHANGELOG.md`, commit, `git tag vX.Y.Z && git push
origin vX.Y.Z`; the release workflow builds and checks both packages and the Windows app and publishes the GitHub
Release; then merge to `main`, which redeploys the website. Details: [docs/development.md](docs/development.md#releasing-coal).
