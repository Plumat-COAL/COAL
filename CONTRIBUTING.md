# Contributing to COAL

COAL is two engines that must behave identically (Python in `src/coal`, C++ in `cpp/`), a `pal2nal` drop-in, and an
app in `apps/coal-app`. This page is the short version of how to work on it; the details are in
[docs/development.md](docs/development.md) and [tests/README.md](tests/README.md).

## Licence

COAL is under GPL-3.0-only ([LICENSE](LICENSE), [NOTICE](NOTICE); 0.5.0 and earlier were MIT). By contributing you
agree that your contribution is distributed under GPL-3.0-only. Code copied or adapted from elsewhere must be under
a licence compatible with GPL version 3 and must keep its copyright notice; say where it comes from in the file and
add it to `NOTICE`. PAL2NAL's code (GPL-2.0) must not be copied in.

## Set up

```
git clone https://github.com/Plumat-COAL/COAL && cd COAL
python -m venv .venv && . .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -e . -e "apps/coal-app[gui]"      # the engine (coal-align) and the app (coal-align-app)
python -m pip install psutil                            # tests only (the worker-process tests)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOAL_TEST_HOOKS=ON && cmake --build build --config Release -j 4   # coal-cpp, a test build (docs/development.md "Test hooks")
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
| Engine, C++ engine included: every module in `tests/` with the binary (goldens with `COAL_IMPL=cpp`, C++ vs Python, CLI, readers, tables, thread pool, Ka/Ks, alignment layer, pipeline, and every other test that skips without `coal-cpp`), failing on a "not built" skip and on any skip in the core C++ modules but the platform ones (CI's C++ job, on all three systems; `-DCOAL_TEST_HOOKS=ON` build, DIAMOND installed) | `COAL_IMPL=cpp COAL_CPP_BIN=$PWD/build/coal-cpp COAL_REQUIRE_CPP=1 COAL_REQUIRE_ALIGN=1 python .github/ci_cpp_tests.py --no-skip test_differential,test_cpp_vs_python,test_cpp_cli,test_cpp_readers,test_cpp_tables,test_cpp_pool,test_wheels --allow "<CI's matrix.allow for your system: empty on Linux>" discover -s tests` (on Linux CI adds `test_kaks_dispatch` to `--no-skip`) |
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

## Nightly checks

`.github/workflows/nightly.yml` runs every night at 03:17 UTC (and from Actions > Nightly > Run workflow, with an
optional `seed`) the random campaign that found real bugs, with a new seed each night (the run number; the seed is
in the run's title), and the tutorials on their real data. Seven jobs, each one of the commands below, so a failing
night can be replayed on a laptop with the same seed:

| Job | What | Locally (after `python scripts/fetch_pal2nal.py` and `cmake -S . -B build && cmake --build build`) |
|---|---|---|
| `pal2nal-fuzz` | 3,000 random cases x 4 formats against the real PAL2NAL v14, once with the Python engine and once with `coal-cpp` | `PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl FUZZ_N=3000 FUZZ_SEED=<seed> FUZZ_KEEP=fuzz-keep python tests/fuzz_vs_pal2nal.py -v`, then the same with `COAL_FUZZ_ENGINE=cpp COAL_CPP_BIN=build/coal-cpp` |
| `engine-fuzz` (Linux, macOS) | `coal-cpp` vs the Python engine on 1,500 cases of `tests/tools/engine_fuzz.py`'s wider input space (malformed files, every option, reports, `-o`, `--json`), bytes compared | `COAL_CPP_BIN=build/coal-cpp python tests/tools/engine_fuzz.py --n 1500 --seed <seed> --keep engine-fuzz-keep` |
| `kaks-properties` | the Ka/Ks engine's property fuzz (`tests/test_kaks_properties.py`): the unit tests with `COAL_SLOW=1` (1,500 pairs) and a 1,500-pair campaign (properties, swap / permutation invariance, both formats, the Python command, determinism); a finding that is not a known bug (`KNOWN`) fails, and the shrunk cases are uploaded (`kaks-properties-seed-<seed>`) | `COAL_CPP_BIN=build/coal-cpp COAL_SLOW=1 KAKS_FUZZ_SEED=<seed> python -m unittest tests.test_kaks_properties` and `python tests/test_kaks_properties.py --n 1500 --seed <seed> --keep kaks-props-keep` |
| `sanitizers` | `coal-cpp` built with ASan (leaks included) and UBSan, warnings as errors: the C++ job's test set (every module in `tests/` but `test_wheels`; with DIAMOND, Clustal Omega, MAFFT and psutil installed), then 300 generator cases and 300 Ka/Ks property pairs (`python tests/test_kaks_properties.py --cpp build-asan/coal-cpp --n 300 --seed <seed> --max-codons 1000`); any sanitizer report fails | `cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DCOAL_WARNINGS_AS_ERRORS=ON -DCOAL_TEST_HOOKS=ON -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer -O1 -g" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" && cmake --build build-asan`, then `ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 COAL_SKIP_MEMORY_TESTS=1 COAL_IMPL=cpp COAL_CPP_BIN=$PWD/build-asan/coal-cpp python .github/ci_cpp_tests.py --exclude test_wheels discover -s tests` (the C++ job's test set, every module in `tests/` but the wheel test) and `python tests/tools/engine_fuzz.py --n 300 --seed <seed> --sanitizer` with the same `COAL_CPP_BIN` |
| `tsan` | `coal-cpp` and `coal-cpp-align` built with ThreadSanitizer: `coal-cpp kaks`, `dnds` and `codon-dnds -j 4` with all nine methods (Comeron and Li included) in the default and the reference flavours, and the aligner with `-j 4`; a ThreadSanitizer report fails, and so does an output that differs from the same command's with `-j 1` | `cmake -S . -B build-tsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCOAL_WARNINGS_AS_ERRORS=ON -DCMAKE_CXX_FLAGS="-fsanitize=thread -g -O1 -Wno-error=free-nonheap-object" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" && cmake --build build-tsan`, then `python tests/tools/tsan_threads.py --build build-tsan --seed <seed>` (about 15 minutes on two cores; if ThreadSanitizer stops with "unexpected memory mapping", `sudo sysctl vm.mmap_rnd_bits=28`) |
| `cpp-extras` | the coal-cpp tests the C++ job of ci.yml skips for want of something other than the binary: `test_cpp_vs_reference` (the fuzz against PAL2NAL v14, the hookless-build test), `test_stress_findings` (`COAL_SLOW`), `test_align` (Biopython), `test_robustness.TestDiskFull` (a full tmpfs, as root); none may skip | `COAL_IMPL=cpp COAL_CPP_BIN=$PWD/build/coal-cpp COAL_CPP_NOHOOKS=$PWD/build-nohooks/coal-cpp COAL_SLOW=1 PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python .github/ci_cpp_tests.py --no-skip test_cpp_vs_reference,test_stress_findings,test_align tests.test_cpp_vs_reference tests.test_stress_findings tests.test_align` (`build-nohooks`: a build without `-DCOAL_TEST_HOOKS=ON`; `pip install biopython`), and as root `python .github/ci_cpp_tests.py --no-skip test_robustness tests.test_robustness.TestDiskFull` |
| `tutorials` | ci.yml's tutorials job without its paths filter: the commands of `docs/tutorials/` 01, 02, 03 and 05 on their real NCBI data (with tutorial 3's MAFFT run), compared with the outputs, exit statuses and numbers the tutorials show; catches drift from a new tool pin or from NCBI's files. Not seeded | `python scripts/run_tutorials.py --full --jobs 4` (needs `COAL_CPP`, DIAMOND, Clustal Omega and MAFFT: `coal pipeline install diamond clustalo mafft`; PAL2NAL v14 optional) |

`COAL_SKIP_MEMORY_TESTS=1` skips the peak-memory test (`tests/test_cpp_vs_python.py` `TestCppMemory`), which is
meaningless under ASan's shadow memory. The engine generator has its own quick test, `tests/test_engine_fuzz.py`
(20 cases; skips without `coal-cpp`), so the tool cannot rot between nights.

**Reading a failure.** A failed fuzz job uploads an artifact (`pal2nal-fuzz-seed-<seed>`,
`engine-fuzz-<os>-seed-<seed>`, `kaks-properties-seed-<seed>`, `sanitizers-seed-<seed>`, `tsan-seed-<seed>`,
`tutorials-nightly`; kept 30 days) with the step's full log (`fuzz-*.txt`,
`engine-fuzz.txt`, `san-*.txt`: the `DIFF #n` / `SANITIZER #n` lines say what differed and show both outputs) and
one folder per differing case: `fuzz-keep/python/c<n>/` or `fuzz-keep/cpp/c<n>/` (PAL2NAL fuzz; `ref*/` inside are
the neutralised reference copies), `engine-fuzz-keep/c<n>/` (generator), `san-keep/san<n>/` and `san-kaks-keep/` (sanitizer hit),
`kaks-props-keep/` (Ka/Ks property finding), and the tutorials' `report.json` and logs. Each
folder holds the input files and `args.txt`, the command line that was run: `cd` into it and run
`pal2nal $(cat args.txt)` / `perl .../pal2nal.pl $(cat args.txt)` (PAL2NAL fuzz) or `coal $(cat args.txt)` and
`coal-cpp $(cat args.txt)` (generator) to reproduce. An unexplained difference is a bug in COAL (or a missing entry
in `DELIBERATE` / `docs/differences.md`) until shown otherwise; a sanitizer report is a bug in `cpp/src`.

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

The release checklist is in one place, [docs/development.md, "Releasing COAL"](docs/development.md#releasing-coal):
the version in the six places `scripts/release_version.py` checks (and the others it lists), the CHANGELOG
section, the generated files, the coalr install lines, the prerelease differences, the tag, and the next
development version afterwards. `python scripts/release_version.py vX.Y.Z` checks a commit against a tag the way
`.github/workflows/release.yml` does (exit 1 with the mismatch spelled out); `--dry-run` only checks that the six
places agree.
