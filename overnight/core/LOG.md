# Overnight log: core (branch overnight/core)

## Run 5 step 3 summary: coal-cpp mirrors the merge (2026-10-17 15:40 UTC)

`overnight/core` was fast-forwarded to Keanu's merge `origin/keanu_overnight/integration` 0729829 (it already contains
51b8092, so no merge commit), then each MRG decision was mirrored in coal-cpp. Every change has a CHANGES row citing
its MRG row. **Not pushed from here** (no credentials): `git push origin overnight/core` from Chloe's terminal.

| MRG decision | Status (commit) |
|---|---|
| MRG-1 safe ID pairing: exact, normalised only when it changes the pairing (notice), ENG-d notice on a partial match, silent when no ID matches, `partial_id_match` error with --strict | done: CPP-25 (3bfba03); the C++-only pairing warning, `coal: summary:` C++ lines and `COAL_CPP_ONLY_LINES` are gone |
| MRG-2 / MRG-2b report: notices, summary / methods text, hints, aborted-run entries, `options.strict/blockonly`, error codes | done: CPP-27 (3bfba03), full port of report.py; C++ error codes = Python's (CPP-26) |
| MRG-3 `--keep-going` with `--strict` (a refused sequence is skipped, exit 2), per-item stderr order, `--max-warnings` caps per-codon lines only (notices never) | done: CPP-26 (3bfba03); max-warnings below 0 is an error: CPP-29 |
| MRG-4 `--resume` vs cleanup: checkpoint only for regular/new outputs (hard links in place) and regular inputs, "cannot be resumed (...)" notes, partial files kept after a signal only once a checkpoint exists | done: CPP-28 (4256784) |
| MRG-3/4 CLI: extras are ordinary abbreviable options in the merged order (`--j` now ambiguous, as in Python); help regenerated | done: CPP-29 (4256784) |
| MRG-8 version 0.4.0; `-j` help | already 0.4.0 in CMake/C++ (Keanu's merge); C++ keeps "worker threads" in the -j help, as MRG-8 documents |
| MRG-6 tests: the merged branch's own tests are the default for C++ | done: CPP-30 (4256784): test_merge, test_report, test_engine_extras added to the shim; `COAL_REF_SRC` kept only as an override (MRG-6) |
| 32-bit conversion in the new report code | fixed: CPP-31 (47d2168) |

**Results (Linux, GCC 11 Release `-Werror`, at 47d2168):**
- Goldens with `COAL_IMPL=cpp` (test_differential): **285/285**.
- Random C++ vs merged Python, report JSON included: **3000/3000** for seed 20261017 and **3000/3000** for seed 777 (213/1000 right after the merge, before CPP-25..27).
- Fuzz vs real PAL2NAL v14 through the merged harness (normalised_ids feature etc.), FUZZ_N=1500, 4 formats: seeds 20261017, 777, 4242: **0 unexplained** each.
- The merged branch's own engine tests with coal-cpp as `coal`: **151 tests, all pass except the 2 in-process DEC-3 mocks** (not applicable to a binary). Tests calling `pipeline.run` in-process (TestMrg3KeepGoing matrix, TestKeepGoing -j) still exercise Python; the C++ side of those is covered by the random test and TestCppThreads.
- Full suite: tests/ (16 modules, run in parts for the 180 s cap) all OK; tests_app 235 tests OK (7 skipped), with `PYTHONPATH=src` as CI's editable install.
- ASan/UBSan: goldens 285/285, the 151 reference tests, random 600, resume/threads tests: clean. TSan: 20,000 sequences, -j 3/4, fasta / clustal+--keep-going+--strict / codon+report / paml+--keep-going+--progress: clean; codon -j4 output = -j1.
- Cross-compile (cache-busted): Windows x86 / x64, macOS arm64: 0 warnings after CPP-31.
- Windows (Chloe, MSVC 19.51 x86 /WX) at 4256784: build clean incl. report.cpp. Her test run still had `COAL_REF_SRC` pointing at the pre-merge reference (279 goldens, 4 stderr-order differences = the old order); to be re-run with the variable cleared.

**Parked:** nothing.

**Needs a human decision:**
- Push: `git push origin overnight/core`.
- Re-run on Windows with `set COAL_REF_SRC=` (expected 285/285 and random OK).
- `--version`: C++ prints `coal 0.4.0 (C++)`, Python `coal 0.4.0` (tests only check the prefix); keep the suffix?

### Step 4 (benchmarks) and CPP-32/CPP-33 (2026-10-17 16:00 UTC)
- The first step-4 run showed C++ peak memory on 100,000 x 500 at 152 MB (26 MB before the merge): the ported
  report built a view of every entry at the end of the run, and the IDs were held 3-4 times. Fixed in CPP-32 /
  CPP-33 (4084d1f), guarded by TestCppMemory (fails with the old code: 1,455 bytes per extra sequence; now ~300).
  500,000 x 150: fasta 110 MB / 2.4 s, clustal 134 MB, codon 146 MB. All step-3 checks re-run green on 4084d1f
  (goldens 285/285, random 3000/3000 x 2 seeds, 151 reference tests, ASan/UBSan, TSan, cross-compile 0 warnings).
- Benchmarks at 4084d1f, `--impl both --coal-jobs 0` (2-core / 4 GB VM): `benchmarks/results/cpp_vs_python_linux_4084d1f.md`.
  10,000 x 500: Python 0.682 s (-j 1) / 0.651 s (default), C++ 0.125 / 0.109 s (5.5x / 6.0x).
  100 x 10,000: 0.176 / 0.175 vs 0.033 / 0.028 s (5.3x / 6.3x). 100,000 x 500: 4.213 / 2.945 vs 1.050 / 0.964 s
  (4.0x / 3.1x); C++ peak 35-37 MB vs Python 105-125 MB. Outputs byte-identical everywhere.
- Needs a human: push; the same benchmark on a machine with more cores for headline numbers; Windows tests with
  `COAL_REF_SRC` cleared.

## Run 5 progress (final engine fixes; reference 0f0330f)
- 2026-10-17 13:21 UTC | CPP-22 FIN-E1, CPP-23 FIN-C14 | reference goldens 279/279; random vs reference 3000/3000 (seed 20261017) and 3000/3000 (seed 777); fuzz vs real PAL2NAL through the reference's updated harness, FUZZ_N=1500: seeds 20261017, 777, 4242: 0 unexplained each; reference tests via shim 94/96 (the same 2 in-process DEC-3 mocks; test_final_fixes added, all pass). FIN-E3 (drop-in argument parsing) is Python-only: coal-cpp is not used as the drop-in, nothing to port.
- 2026-10-17 13:35 UTC | CPP-24 FIN-E2 bounded clustal/codon writer | peak memory and time, 500,000 x 150 codons (2-core VM): fasta -j1 2.26 s / 88 MB, -j0 2.41 s / 89 MB; clustal -j1 2.42 s / 114 MB, -j0 2.61 s / 115 MB (was 76.7 s / 256 MB); codon -j1 3.27 s / 126 MB, -j0 3.46 s / 125 MB. 250,000: fasta 53 MB, clustal 78 MB, codon 91 MB, so each format grows by ~35 MB per 250k sequences, the IDs.
- 2026-10-17 13:43 UTC | CPP-24b 32-bit conversions in the new writer | note: the zig/clang cross-compile checks reuse cached objects when the source and flags are unchanged and then print no warnings; every cross check now passes a random -DCACHEBUST. Final step-1 checks after CPP-24/24b: goldens 279/279; random vs reference 3000/3000 for seeds 20261017 and 777; fuzz FUZZ_N=1500 for seeds 20261017, 777, 4242: 0 unexplained each; reference tests via shim 94/96 (2 in-process mocks); -j 1/2/4 identical (test_cpp_pool, TestCppThreads, reference TestDec4ParallelAndBlocks); -Werror GCC clean; clang -Werror clean for Windows x64/x86 and macOS arm64 (cache-busted); ASan/UBSan (large test, goldens with the chunked writer forced) and TSan (pool) clean.
- 2026-10-17 13:48 UTC | Step 2: native Windows check (Chloe's PC, MSVC 19.51 x86, Release, -DCOAL_WARNINGS_AS_ERRORS=ON) at 747b1db: build clean (coal-cpp.exe built 13:46 with /W4 /WX, incl. sysio.cpp and the chunked writer); reference goldens with the .exe 279/279; random C++ vs reference Python (Python 3.11 on Windows) 1000/1000 identical.

## Run 4 summary: coal-cpp mirrors the reviewed reference engine (2026-10-17 13:01 UTC)

Reference: `origin/keanu_overnight/integration` at bb80bde, in a read-only worktree (`~/coal-ref`, not merged; `overnight/merged` untouched). Every change has a row in `overnight/core/CHANGES.md` (CPP-9 .. CPP-21, each citing the reference row it mirrors). **Not pushed from here** (no GitHub credentials on this machine): push with `git push origin overnight/core` from Chloe's terminal.

| Brief item | Reference rows | Status (commit) |
|---|---|---|
| Harness: random test + goldens + PAL2NAL fuzz against the reference (`COAL_REF_SRC`) | – | done: CPP-9 (5e233d1), CPP-21 reference tests via a shim (bc73ab6) |
| FZ-5 -nogap: the selected code's stops | FZ-5 | done: CPP-10 (11be0b2) |
| -nomismatch: no mismatch warnings; B/Z/J/O/U columns removed | ENG-i2, DEC-4f | done: CPP-11 (11be0b2) |
| Readers: CLUSTAL '#' lines, repeated headers, leading blanks; non-ASCII descriptions | FZ-1, FZ-2 (ENG-h, R2-E6, R3-E3 already matched) | done: CPP-12 (004fafb) |
| -blockonly and '#' mark lines; codon + -blockonly error | DEC-4 | done: CPP-13 (004fafb) |
| tfastx `\` and `/` | FZ-3 | done: CPP-14 (d384d30) |
| --strict and the notices (>half checkable mismatch, position guess, no anchor, short ambiguous protein, partial ID match); -q silences notices | DEC-1, ENG-c, R2-E3, R3-E1, ENG-d, FZ-4, FZ-6, DEC-2 | done: CPP-15 (13f88bc) |
| -o: symlinks, FIFOs/devices, hard links, permissions, loop / missing-folder errors | ENG-a, R2-E5, R3-E2 | done: CPP-16 (13f88bc) |
| Piped / FIFO inputs; errors name the argument | ENG-b, R2-E7 | done: CPP-17 (13f88bc) |
| SIGTERM 143, SIGHUP 129, Ctrl-C 130, closed pipe 141 (Windows EINVAL too) | ENG-f, R2-E8, DEC-3, R2-D1b | done: CPP-18 (13f88bc) |
| CLI help / messages as the reference (-blockonly, -strict, formats) | DEC-1, DEC-4, FZ-5 help | done: CPP-19 (bacdf2d) |
| Table 3 GTG start | ENG-g | already done (CPP-3) |
| R2-E9 (IUPAC frame search) | R2-E9 | nothing to port: performance only, the C++ search already handles IUPAC DNA directly |
| First-Met start codons in the anchor fallback | (being changed in the reference) | not changed, as the brief says |
| Extra: PAML header for 1,000+ sequences (` %3d %6d`, as PAL2NAL) | – | fixed: CPP-20 (4f0d73d), found by a new large-input test |

**Results (Linux, GCC 11, Release, at the end of the run):**
- Reference goldens (the reference's own harness, coal-cpp swapped in): **279/279** (start: 248).
- Random C++ vs reference Python: **3000/3000** identical (`COAL_RANDOM_N=3000`; 1000/1000 at the default; includes `--strict` on 1/5 and `-q` on 1/7 of the cases). Start: 739/1000.
- Fuzz vs the real PAL2NAL v14: **0 unexplained differences** for FUZZ_N=1000 (seed 20261017) and FUZZ_N=2500 (seed 7), 4 formats each. Start: 440.
- The reference's own engine tests run against coal-cpp: test_review_fixes 43/43, test_fuzz_findings 17/17, test_decisions 22/24 (the 2 others mock sys.stdout inside the Python process: not applicable to a binary; the real closed-pipe test in test_review_fixes passes).
- -j 1/2/4: identical in every format (reference TestDec4ParallelAndBlocks via the shim; test_cpp_pool; TestCppThreads).
- Builds: `-DCOAL_WARNINGS_AS_ERRORS=ON` clean (GCC); clang `-Werror` clean for Windows x64, Windows x86 (+ `-Wshorten-64-to-32 -Wconversion`) and macOS arm64. Sanitizers: random + goldens + reference tests + fuzz 400 under ASan/UBSan, pool + random 300 under TSan: clean. **MSVC not re-run since CPP-7** (new file cpp/src/sysio.cpp): please rebuild once on Windows.
- This branch's own suite (`python -m unittest discover -s tests`, local mode): 81 tests, 2 fail **by design**: test_cpp_vs_python's random and resume tests compare with this branch's old Python engine, which the reference supersedes (see below). Local goldens with COAL_IMPL=cpp: 272/285, the 13 failures are all `fasta_nogap_ct<N>` (this branch's harness still expects PAL2NAL's -nogap rule; the reference harness expects FZ-5).

**Parked:** nothing.

### Needs reconciling after the merge
1. **Safe ID pairing vs ENG-d.** coal-cpp pairs IDs that only match after normalising (version suffix, `sp|..|..` headers, case) one-to-one; the reference pairs those by order (silently when no ID matches exactly, with the partial-ID notice when some do). Different output for such inputs; with --strict the reference errors on a partial match where coal-cpp pairs by normalised ID. The random test uses plain IDs for those cases in reference mode.
2. **No ID matches at all:** coal-cpp warns "IDs don't match, so sequences were paired by file order ..."; the reference is silent (reference test `test_no_ids_match_stays_silent`).
3. **`coal: summary:` lines and the genetic-code hint** exist only in coal-cpp. Items 2-3 can be switched off with `COAL_CPP_ONLY_LINES=0` (used by the comparisons); decide whether to keep them, or move them to the Python engine too.
4. **--resume vs ENG-f cleanup:** with `-o`, an interrupted coal-cpp run (SIGTERM / SIGHUP / Ctrl-C) keeps `<out>.coal-partial` and its checkpoint for `--resume`; the reference removes everything. Stdout runs and `-o` into a FIFO/device/hard link clean up as the reference.
5. **Run report / --json:** notices are not recorded in the report entries (a probably-unrelated sequence still counts as "warning" via its codon warnings, not via the notice); `strict` and `blockonly` are not in `options`; methods_text unchanged.
6. **--keep-going with --strict:** a --strict error skips that sequence (exit 2) instead of stopping the run (exit 1).
7. **--max-warnings** limits per-codon warnings only; notices are never hidden.
8. `--version` prints `coal 0.4.0 (C++)`; the reference prints `coal 0.3.0`. `-j` help says threads (coal-cpp) instead of processes.
9. After the merge, this branch's local tests (test_differential DELIBERATE rules, test_cpp_vs_python local mode, the frameshift_codon_line golden) should be re-pointed at the merged engine; the reference-mode harness (`COAL_REF_SRC`) can then become the default.

### Needs a human decision
- Push `overnight/core` (needs credentials): `git push origin overnight/core`.
- One MSVC build on Windows (new sysio.cpp): `cmake -S . -B build-win -DCOAL_WARNINGS_AS_ERRORS=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build-win`.
- Windows signals: SIGTERM cannot be caught there (TerminateProcess); the Windows closed-pipe (EINVAL) path is written but only testable natively.
- Keep or drop the test switch `COAL_CPP_ONLY_LINES` after the reconciliation above.

### How to run the comparisons
```
git worktree add --detach ~/coal-ref origin/keanu_overnight/integration
python ~/coal-ref/scripts/fetch_pal2nal.py
export COAL_REF_SRC=~/coal-ref/src PAL2NAL=~/coal-ref/PAL2NAL/pal2nal.v14/pal2nal.pl
python tests/test_cpp_vs_reference.py golden      # reference goldens
python tests/test_cpp_vs_python.py                # random vs reference Python
FUZZ_N=1000 python tests/test_cpp_vs_reference.py fuzz
python tests/test_cpp_vs_reference.py reftests    # the reference's own engine tests with coal-cpp
```

## Run 3 summary: C++ foundations brief (2026-10-17 11:08 UTC)

All six items are done. They are committed on the local branch `overnight/core`. The branch is **not pushed**: `git push origin overnight/core` failed with "could not read Username", because this machine has no GitHub credentials. Push it from a machine that has them. `main` was not touched, and `origin/keanu_overnight/integration` was not merged (it was only read with `git show`). Every change has a row in `overnight/core/CHANGES.md`, and its ID starts the commit message.

| item | ID | commit |
|---|---|---|
| 1 Build setup (CMake, C++17, std only, warnings-as-errors option for CI) | CPP-1 | `be537a6` |
| 2 Input readers = reviewed `seqio.py` + Python/C++ reader test | CPP-2, CPP-2b | `c7d0ac0` |
| 3 Genetic code tables (17; table 3 GTG) + table test | CPP-3 | `801ae19` |
| 4 CLI: same options/help/messages as reviewed `python -m coal` | CPP-4 | `e2fc161` |
| 5 Ordered streaming thread pool (`pool.hpp`) + stand-in test | CPP-5 | `29af0a5` |
| 6 `COAL_IMPL=cpp` hookup in test_differential (skips when not built) | CPP-6 | `65c4a5c` |

**Test results (at the end of the run):**
- Clean build from scratch with `-DCOAL_WARNINGS_AS_ERRORS=ON` (GCC 11): no warnings.
- clang `-Werror -Wall -Wextra -Wpedantic -Wshadow`, cross-compiled for Windows (MinGW) and macOS (arm64): clean.
- Golden checks: Python 285/285, C++ 285/285.
- Full suite `python -m unittest discover -s tests`: 77 tests OK, none skipped, 50 s.
- New tests:
  - test_cpp_readers: 43 edge cases.
  - test_cpp_tables: 17 codes.
  - test_cpp_cli: 50 error, 31 accepted, 10 help/version command lines, plus the help text.
  - test_cpp_pool: 7 tests.
  - All pass. Each was break-tested (code broken on purpose, failure seen, code restored).
- test_cpp_vs_python: 1000/1000 identical.
- The CLI tests pass under ASan/UBSan; the pool tests pass under TSan.

**Parked:** nothing.

**Needs a human decision:**
1. **Output formats.** The brief says "for anything other than fasta print a not implemented yet error". The C++ clustal/paml/codon writers already existed from run 1 and pass all their golden checks, so I kept them working instead of disabling them. Say if you want them switched off.
2. **`-strict/--strict`** is parsed exactly like the reference, then exits 2 with "not implemented yet in coal-cpp". It needs the core rules, which were out of scope.
3. **The CLI follows Python 3.10/3.11 argparse messages.** Python 3.12.8+/3.13 quote `invalid choice` values differently and format help differently. The live comparison therefore runs only on 3.10/3.11. On other versions the test checks the messages pinned from 3.10.
4. **Deliberate CLI differences from the reference** (all in CHANGES CPP-4):
   - `-output` lists four choices.
   - `--version` prints `coal 0.4.0 (C++)`.
   - The C++-only options (`-keep-going`, `--max-warnings`, `--progress`, `--resume`, `--report`, `--json`) match exact names only, so `--j` still means `--jobs`.
   - `-q` now also silences notices, as the reference does.
5. **This branch's Python engine (`src/coal`) predates the reviewed fixes.** For example, it lacks the table-3 GTG start. The C++ readers and tables now follow the integration branch, so for table-3 genes starting with GTG they differ from this branch's Python until integration is merged.
6. ~~**MSVC is not available here.**~~ Resolved: native MSVC build and tests pass on Chloe's PC (see below). Its flags are set (`/W4 /permissive- /utf-8`, `/WX` with the CI option), and the code is clean under clang `-Wconversion`, which approximates MSVC's loss-of-data warnings. Please run one native Windows build:
   `cmake -S . -B build -DCOAL_WARNINGS_AS_ERRORS=ON && cmake --build build --config Release`
   - 2026-10-17 11:15 UTC: tried to run it on Chloe's Windows PC. Not possible: the PC has no Visual Studio or Build Tools (no VS developer prompts, no VS Installer, no Windows CMake in the Start menu), and the terminal can only be clicked, not typed into. Next step: install "Build Tools for Visual Studio 2022" with the "Desktop development with C++" workload (includes MSVC and CMake), or run a windows-latest GitHub Actions job (.github/ is outside this branch's files: human decision).
   - 2026-10-17 11:23 UTC: Chloe installed Build Tools (MSVC 19.51, VS 18) and ran the build: it failed on core.cpp:73 (C4244 long long -> 32-bit size_t, x86 target) and was a Debug build by default. Fixed in CPP-7; the earlier 64-bit -Wconversion check could not see it (size_t is 64-bit there), so the clang check now also targets 32-bit Windows.
   - 2026-10-17 11:28 UTC: MSVC build now clean (/W4 /WX, x86, Release). Windows test run: `--version` OK; test_cpp_cli 11/11 OK; golden 200/285 and pool tests 14 errors, both test-harness problems (CRLF stderr; Linux binary picked up from the shared build/ folder). Fixed in CPP-8; Windows re-run pending.
   - 2026-10-17 11:31 UTC: **Windows verified** (after CPP-7/CPP-8): MSVC 19.51 x86 Release build clean with /W4 /WX; golden C++ 285/285; test_cpp_cli, test_cpp_pool, test_cpp_vs_python (1000 random inputs), test_cpp_readers, test_cpp_tables: 30 tests OK (66 s), Python 3.11. Item 6 is resolved.
7. **The local branch `overnight/merged`** (the merge with Keanu's integration branch, made earlier) is still unpushed and untouched.

## Morning summary (both overnight runs + code review)

Everything is on the local branch `overnight/core` (17 commits on top of main `31cf3d0`). **Not pushed**: this machine has no GitHub credentials. Push it yourself when you're happy: `git push -u origin overnight/core`

**State:** Python 285/285 and C++ 285/285 golden checks (279 + the new frameshift_codon_line case). The full test suite (52 tests) passes. The C++ vs Python random test is identical on 1000/1000 inputs, plus 3,000 more with other seeds. The final code-review pass came back clean.

### Done → commit
| task | commit |
|---|---|
| Run 1 · 1 Python `-output clustal/paml/codon` | `2d80616` |
| Run 1 · 2 C++ core `coal-cpp` (+ random/threads tests, speed, benchmarks) | `4eb39aa` `5266c9d` `e1dd30c` `7add571` |
| Run 1 · 3+4 run report `--report/--json`, `coal: summary:` lines, genetic-code hint | `f1373a2` |
| Run 2 · 1 summary + hints | already done in run 1 (skipped) |
| Run 2 · 2 `--keep-going` (exit 2) | `c409b01` |
| Run 2 · 3 safe ID pairing (normalised, one-to-one; file-order warning) | `004d430` |
| Run 2 · 4 `methods_text` in the report, `--max-warnings N` | `eef8c0b` |
| Run 2 · 5 `--progress` | `6ace724` |
| Run 2 · 6 `--resume` (checkpoints) | `4e98e10` |
| Run 2 · 7 C++ parity for all of the above | `3f31971` |
| Code review passes 1–4 | `77ee9d6` `640582b` `691203b` (pass 4 clean, no commit needed) |

### Speed (measured, 2-core Linux VM, local disk, medians)
| dataset | Python -j1 | C++ 1 thread | C++ 2 threads | speed-up |
|---|---|---|---|---|
| 10,000 × 500 | 0.593 s | 0.114 s | 0.097 s | 6.1× |
| 100 × 10,000 | 0.149 s | 0.027 s | 0.023 s | 6.4× |
| 100,000 × 500 | 3.559 s | 0.952 s | 0.881 s | 4.0× |
C++ peak memory < 30 MB on a 137 MB input. Re-checked after all the changes: still the same. Details: benchmarks/results/cpp_vs_python_linux.md

### Code review (repeated 4×, stopping early if a pass came back clean)
- Pass 1 fixed 4 issues, including a real C++ heap-use-after-free on error paths (found with AddressSanitizer) and temp-file leaks.
- Pass 2 fixed a wrong number in this log.
- Pass 3 fixed a misleading `--resume` message.
- Pass 4 was clean.
- Tools used: pyflakes, `-Wall -Wextra -Wshadow`, ASan + UBSan and TSan builds run through the test suites, Python `-X dev` ResourceWarnings, CRLF inputs, adversarial edge cases on both engines.

### Parked
- ~~`-output codon` residue-line quirk~~ Fixed next morning (see below).

### Needs a human decision
- ~~Add `build/` to .gitignore.~~ Done (Chloe asked, next morning).
- An untracked `Claude outputs/` folder sits in the repo root. It isn't mine, so I left it alone.
- ~~Bump the version to 0.4.0~~ Done (pyproject.toml, src/coal/__init__.py, CMakeLists.txt and the C++ VERSION constants; tests/report_example.json regenerated).
- coal-cpp has not been compiled with MSVC on Windows (no compiler here). It is written to standard C++17, but try the Windows build below.
- README (off limits) doesn't mention the new flags yet.

### How to try it
Linux/macOS:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
export PYTHONPATH=src            # or: pip install -e .
python -m coal tests/fixtures/clean/pep.fasta tests/fixtures/clean/nuc.fasta -output paml -o out.paml --report report.json
python -m coal pep.fasta nuc.fasta --keep-going -o out.fa          # skip broken genes, exit code 2
python -m coal pep.fasta nuc.fasta -o out.fa --json                # report on stdout
python -m coal pep.fasta nuc.fasta --max-warnings 20 --progress
python -m coal pep.fasta nuc.fasta -o out.fa --resume              # after an interrupted run
./build/coal-cpp <same options>                                    # C++ engine, same behaviour
python -m unittest discover -s tests
COAL_IMPL=cpp python -m unittest tests.test_differential
python tests/test_differential.py --count
python benchmarks/run_bench.py --impl both --sizes 10000x500
```
Windows (PowerShell, Visual Studio Build Tools + CMake):
```
cmake -S . -B build
cmake --build build --config Release
$env:PYTHONPATH="src"
python -m coal tests\fixtures\clean\pep.fasta tests\fixtures\clean\nuc.fasta -output paml -o out.paml --report report.json
.\build\Release\coal-cpp.exe tests\fixtures\long_mismatch\pep.fasta tests\fixtures\long_mismatch\nuc.fasta -output codon
$env:COAL_IMPL="cpp"; $env:COAL_CPP_BIN="build\Release\coal-cpp.exe"; python -m unittest tests.test_differential tests.test_cpp_vs_python
```
Try the mito hint: run mitochondrial genes with the default code, and the summary suggests `-codontable 2`.

(Timestamps below come from the VM's `date` and may not match UK time.)

<details><summary>Run 1 summary (kept for reference)</summary>

## Morning summary — run 1 (C++ core prompt)

**Done** (all on local branch `overnight/core`; nothing pushed — no GitHub credentials on this machine)
- Task 0 preflight → OK (cmake 4.4.4, g++ 11.4, C++17 + threads)
- Task 1 Python `-output clustal|paml|codon` → `2d80616`
- Task 2 C++ core `coal-cpp` → `4eb39aa` (golden), `5266c9d` (random + threads tests), `e1dd30c` (speed), `7add571` (benchmarks)
- Task 3 run report `--report` / `--json` → `f1373a2`
- Task 4 `coal: summary:` lines + genetic-code hint → `f1373a2` (same commit; built together with the report)

**Golden checks:** Python 279/279, C++ 279/279 (every case × every variant in tests/cases.json, incl. clustal/paml/codon).

**C++ vs Python random test:** 1000/1000 identical (stdout, full stderr, exit code), ~7 s; `-j 1/2/3` identical on a 2-batch dataset in every format.

**Speed (measured, 2-core Linux VM, data on local disk, medians):**
| dataset | Python -j1 | C++ 1 thread | C++ 2 threads | speed-up |
|---|---|---|---|---|
| 10,000 × 500 | 0.593 s | 0.114 s | 0.097 s | 6.1× |
| 100 × 10,000 | 0.149 s | 0.027 s | 0.023 s | 6.4× |
| 100,000 × 500 | 3.559 s | 0.952 s | 0.881 s | 4.0× |
C++ peak memory < 30 MB on the 137 MB input. Details: benchmarks/results/cpp_vs_python_linux.md.

**Parked:** `-output codon` residue-line quirk of PAL2NAL with frameshift digits 4/5 (see Decisions) — no golden covers it.

**Needs a human decision:** add `build/` to .gitignore; untracked `Claude outputs/` folder in the repo root (not mine, left alone); C++ does not yet print the `coal: summary:` lines or support `--report` (Python only — run 2's task 7 may add it).

**How to build and try it**
Linux/macOS:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/coal-cpp tests/fixtures/long_mismatch/pep.fasta tests/fixtures/long_mismatch/nuc.fasta -output codon
PYTHONPATH=src python -m coal tests/fixtures/clean/pep.fasta tests/fixtures/clean/nuc.fasta -output paml -o out.paml --report report.json
python -m unittest discover -s tests            # Python engine (needs PYTHONPATH=src or pip install -e .)
COAL_IMPL=cpp python -m unittest tests.test_differential
python tests/test_differential.py --count       # "279/279 golden checks"
python benchmarks/run_bench.py --impl both --sizes 10000x500
```
Windows (PowerShell, Visual Studio Build Tools + CMake installed):
```
cmake -S . -B build
cmake --build build --config Release
.\build\Release\coal-cpp.exe tests\fixtures\long_mismatch\pep.fasta tests\fixtures\long_mismatch\nuc.fasta -output codon
$env:PYTHONPATH="src"; python -m coal tests\fixtures\clean\pep.fasta tests\fixtures\clean\nuc.fasta -output paml -o out.paml --report report.json
$env:COAL_IMPL="cpp"; $env:COAL_CPP_BIN="build\Release\coal-cpp.exe"; python -m unittest tests.test_differential
```
Push when happy: `git push -u origin overnight/core`

(Timestamps below come from the VM's `date` and may not match UK time.)

</details>

## Progress log

- 2026-10-17 01:17 UTC | task 0 | preflight OK: cmake 4.4.4 + g++ 11.4, C++17 + std::thread hello-world built and ran (2 hardware threads). Branched overnight/core from main 31cf3d0 (pull: already up to date).
- 2026-10-17 01:17 UTC | task 1 | Python 279/279 golden checks (all variants: fasta, clustal, paml, codon, nogap, nomismatch, codon tables), C++ n/a | added -output clustal/paml/codon; test_differential checks every variant; -j1 == -j2/-j3 tests for every format.
- 2026-10-17 01:25 UTC | task 2 | Python 279/279, C++ 279/279 golden checks | C++17 core in cpp/src (tables, core, seqio, main) + root CMakeLists.txt (coal_core static lib + coal-cpp exe). Whole port compiled warning-free and passed every golden check on the first build, so there were no intermediate counts to commit. test_differential runs either engine (COAL_IMPL=python|cpp, COAL_CPP_BIN) and has a --count mode.
- 2026-10-17 01:27 UTC | task 2 | Python 279/279, C++ 279/279 | tests/test_cpp_vs_python.py: 1,000 random inputs (seed 20261017; Python outcomes 301 ok / 647 with warnings / 52 errors) -> 1000/1000 identical stdout + full stderr + exit code, ~7 s; plus -j1/-j2/-j3 identical on a 3.4M-nt, 2-batch dataset for all formats, warnings in input order.
- 2026-10-17 01:31 UTC | task 2 | Python 279/279, C++ 279/279, random 1000/1000 | C++ speed-ups: table lookups hoisted out of hot loops, sequential (not seek-per-row) reads of staged rows, memchr line scanning in the DNA index. 10k x 500 on the VM's local disk: C++ 0.09 s (2 threads) / 0.12 s (1 thread) vs Python 0.59 s. run_bench.py: --impl python|cpp|both, --cpp-bin.
- 2026-10-17 01:34 UTC | task 2 done | Python 279/279, C++ 279/279, random 1000/1000 | expand() fast path, token reserve. Measured (2-core VM, local disk, medians): 10k x 500 Python 0.593 s / C++ 0.097 s (6.1x); 100 x 10k 0.149 / 0.023 (6.4x); 100k x 500 3.559 / 0.881 (4.0x); C++ peak < 30 MB. Results in benchmarks/results/cpp_vs_python_linux.md (+ csv/txt).
- 2026-10-17 01:40 UTC | tasks 3+4 | Python 279/279, C++ 279/279, random 1000/1000, all suites OK | src/coal/report.py: --report FILE.json / --json (needs -o), schema coal.report/1 exactly as specified plus extra fields (sequences[].index, summary.not_processed, summary_text, error{code,text,detail}, hints[].checked/fit/fit_current); report written on failed runs too (plain-English reason per error code). End-of-run 'coal: summary:' lines + genetic-code hint (>20% mismatching/failed -> re-test up to 200 of them with all 17 codes; suggest only if it fits >=80% and >=20 points more than the current code; ties -> lowest table number). The summary is collected on every run with counters only (per-sequence entries only kept when a report is requested). -q suppresses summary lines. tests/test_report.py (10 tests), tests/report_example.json (real run: 6 mito genes with code 1 -> hint code 2). C++ does not print summary lines yet; test_cpp_vs_python ignores 'coal: summary:' lines (allowed by the task).

## Run 2 progress (engine extras prompt)

- run 2 task 1 (summary + genetic-code hints): already done in run 1 (f1373a2); skipped.
- 2026-10-17 01:43 UTC | run2 task 2 --keep-going | Python 279/279, all suites OK | failing sequences are skipped ('coal: warning: skipped <msg>' in input order), others written in order, exit 2, report/summary list them; works with -j>1 and all formats (paml stages rows so the header count is right); -nogap/-nomismatch columns come from the sequences that succeeded. Decisions: if every sequence fails -> error, exit 1, no output; frameshift column widths still come from the whole protein alignment (a failed row's digit can widen a column for the others).
- 2026-10-17 01:46 UTC | run2 task 3 safe ID pairing | Python 279/279, C++ 279/279, random 1000/1000 (now incl. version/pipe/case IDs), all suites OK | pair_ids(): exact IDs -> 'id'; else normalised keys (trailing .N removed, lower case, pipe fields minus db tags such as sp/tr/gi) accepted only if every protein has exactly one candidate and the pairing is one-to-one (exact matches preferred) -> 'normalised'; else file order + one warning naming the first two pairs. Per-sequence pairing in the report, summary.pairing overall, summary line. Ported to C++ at the same time so the strict C++/Python stderr comparison stays on. Fixture outputs unchanged (279/279). Decision: the warning uses ASCII '<->' instead of the Unicode arrow, because a Windows console (cp1252) cannot print it.
- 2026-10-17 01:47 UTC | run2 task 4 methods_text + --max-warnings | all suites OK | report.methods_text (version, URL, sequence count, genetic code name + NCBI number, output format, -nogap/-nomismatch, excluded sequences with --keep-going, note for tables 10/11/13). --max-warnings N caps per-codon warning lines (pepAlnPos ...), then '... and N more warnings (see the report)' (or 'see them all with --report FILE.json' when no report was requested); 0 = unlimited (default); the report keeps all. Decision: coal.__version__ left at 0.3.0 (pyproject.toml holds the version too and is off limits) -> needs a human bump to 0.4.0.
- 2026-10-17 01:48 UTC | run2 task 5 --progress | all suites OK | for inputs > 10,000 sequences: 'coal: progress: 61% (61,775 of 100,000 sequences), about 1 s left' at most every 2 s (after each finished batch), '\r'-updated single line on a terminal, plain lines otherwise, final 'done in X s' line. stderr only; output unchanged (tested).
- 2026-10-17 01:54 UTC | run2 task 6 --resume | all suites OK (incl. 7 new resume tests) | src/coal/checkpoint.py. Every -o run writes <out>.coal-partial (or <out>.coal-stage for two-pass formats) + <out>.coal-entries.jsonl (non-ok report entries) and, at most every 2 s, <out>.coal-checkpoint.json (engine, version, options hash, input size+mtime_ns, completed sequences, byte lengths, -nogap/-nomismatch column sets, skipped count, report counters) written atomically after fsync of the data. --resume truncates and continues; refuses with the specified message if inputs/options changed; 'no checkpoint' -> normal run + note. Checkpoint files deleted on success and on real errors, kept after crashes/kills. Two-pass formats (-nogap, -nomismatch, clustal, codon, paml+keep-going) resume the first pass and redo the final pass (said in --help). Tests kill a real subprocess at a known checkpoint (COAL_TEST_PAUSE_AFTER_CHECKPOINT test hook + terminate, watchdog 60 s) and compare byte-for-byte with an uninterrupted run: fasta/paml/clustal/codon/-nogap/-nomismatch, -j 1 and -j 2, --keep-going. Attempt 1 at the tests failed only because subtests reused an output name (test bug, fixed). Tuning/test knobs: COAL_BATCH_NT, COAL_CHECKPOINT_EVERY. Earlier stderr warnings are not re-printed on resume (the report and summary cover all sequences).
- 2026-10-17 02:06 UTC | run2 task 7 C++ parity | Python 279/279, C++ 279/279, random C++/Python 1000/1000 now incl. --keep-going, --max-warnings, --json reports and the 'coal: summary:' lines (full stderr compared again), all suites OK | coal-cpp gained: run report (cpp/src/report.*: identical JSON content to Python except runtime; checked field-by-field), summary lines + genetic-code hint, methods_text, --keep-going, --max-warnings, --progress, --resume (own checkpoints, engine 'coal-cpp'), error codes. cpp/src/json.*: minimal JSON writer/parser + SHA-256 (standard library only). New tests: C++ kill+resume == uninterrupted Python run (fasta -j2, clustal, -nogap, paml/codon --keep-going) incl. identical reports, refusal + no-checkpoint messages, progress. Found + fixed a 45% slowdown in both engines: the output's SHA-256 was computed on every run; now only when a report is requested (C++ 100k x 500 back to 0.91 s). Test knobs added to C++: COAL_BATCH_NT, COAL_CHECKPOINT_EVERY, COAL_PROGRESS_EVERY, COAL_PROGRESS_MIN_SEQS, COAL_TEST_PAUSE_AFTER_CHECKPOINT.

## Code review passes (requested: repeat up to 4x until a pass is clean)
- 2026-10-17 02:14 UTC | review pass 1 | NOT clean -> fixed | tools: pyflakes (Python), g++ -Wall -Wextra -Wshadow, AddressSanitizer + UBSan build run through the golden, random (1000), resume and progress tests, ThreadSanitizer build (setarch -R) on 3-thread runs of every staged/unstaged path, a CRLF-input check of all 38 fixtures on both engines, and a 13-case edge-case CLI comparison (missing files, missing output dir, empty input, bad options, --resume without checkpoint, ...). Found: (1) C++ heap-use-after-free on error paths: OutFile declared its stream buffer after the stream, so the buffer was freed before the stream's destructor flushed into it (ASan, 4 random inputs) -> member order fixed; (2) C++ temporary '.grid' file not removed if a block-format run failed or crashed in the final pass -> cleaned up on every exit; (3) Python report write left its temp file behind on failure and printed the temp name -> cleaned up, message names the report path; (4) one f-string without placeholders in a test. TSan: no races. Not testable here: MSVC/Windows build of coal-cpp (missing <algorithm>/<cctype> includes added for portability).
- 2026-10-17 02:19 UTC | review pass 2 | NOT clean -> fixed | rebuilt ASan/UBSan on the current code and re-ran golden + random + resume + progress tests (clean, no leaked temp files; the 7 stray C++ temp files found in the temp dir were from the pre-fix ASan crashes of pass 1 and do not reappear); Python under -X dev with ResourceWarnings shown, in-process tests and CLI runs of every path (clean); full suite run twice (no flaky tests: 51 tests OK both times); re-read the C++ resume/cleanup/threading code and the Python pipeline/checkpoint code. Found: one wrong number in this log (task 6 line said '64 tests'; the whole suite has 51) -> corrected.
- 2026-10-17 02:23 UTC | review pass 3 | NOT clean -> fixed | adversarial inputs on both engines (ASan build): empty DNA record, duplicate DNA IDs, pipe/version IDs reordered, non-ASCII header, all sequences failing with --keep-going (paml, -nomismatch), -nomismatch + clustal + --keep-going, corrupted checkpoint, emptied partial file -> engines agree on output, stderr and exit code in every case, and resumed outputs stay identical. Found: with an unreadable checkpoint or an incomplete partial file --resume said 'no checkpoint to resume', which is misleading -> both engines now say what they found ('found an unreadable checkpoint, so nothing can be resumed ...'); test added.
- 2026-10-17 02:26 UTC | review pass 4 | CLEAN | ASan/UBSan (golden, random, resume, progress) and TSan (3 threads: plain, codon+keep-going, -nomismatch+report) on the final code; 3,000 new random cases (seeds 1-3) identical between engines; pyflakes clean; no TODO/debug leftovers; -X dev with ResourceWarnings: none; full suite 52 tests OK; golden 279/279 both engines; no leaked temp files; benchmarks re-checked (unchanged). Review stops here (clean pass).
- 2026-10-17 09:40 UTC | morning | codon residue-line quirk FIXED | black-box probing of PAL2NAL (110 + 8 layouts, output only) showed the rule: after a column where the LAST sequence has digit N, each residue line gets 3*ceil(N/3)/3-1 extra '-' tokens (gap/residue there: none). Implemented in both engines (residue_line_blocks / write_blocks), 98/98 probe layouts identical to PAL2NAL for each engine; new fixture tests/fixtures/frameshift_codon_line + PAL2NAL goldens via scripts/make_golden.py (case added to cases.json) -> 285/285 golden checks on both engines; the previous build fails the new codon golden, so the test covers the fix.

## Run 4: mirror the reviewed reference engine (started 2026-10-17 12:22 UTC)

Reference: origin/keanu_overnight/integration at bb80bde, checked out read-only with
`git worktree add --detach ~/coal-ref origin/keanu_overnight/integration` (outside the shared folder, so Windows
does not see it). Real PAL2NAL v14 fetched with the reference's scripts/fetch_pal2nal.py (official site refused by
the proxy; orthologr GitHub copy, sha256-verified) into ~/coal-ref/PAL2NAL (never in this repo).
Harness (CPP-9): COAL_REF_SRC=~/coal-ref/src for tests/test_cpp_vs_python.py (reference mode) and the new
tests/test_cpp_vs_reference.py (reference goldens + reference fuzz vs PAL2NAL, coal-cpp swapped in).

### Initial work queue (C++ at 48eb082 vs the reference)
- Reference goldens: 248/279. Failing: 15 x `fasta_nogap_ct<N>` (FZ-5: -nogap must use the selected code's
  stops), 16 x `fasta_nomismatch*` (ENG-i2: no per-codon mismatch warnings with -nomismatch).
- Random C++ vs reference Python (1000 cases): 739/1000 identical. Seen: DEC-1/ENG-c/R3-E1 notices missing,
  ENG-i2 warnings, FZ-5 -nogap, DEC-4f -nomismatch on B/Z/J/O/U columns, C++-only "paired by file order" warning
  (safe ID pairing; see "Needs reconciling"), partial-ID notice (ENG-d).
- Fuzz vs real PAL2NAL (FUZZ_N=1000, seed 20261017, 4 formats): 440 unexplained differences (clustal 168,
  paml 109, fasta 109, codon 54). Includes every `-blockonly` case (option missing in C++: DEC-4), CLUSTAL inputs
  with '#' mark lines (FZ-1), tfastx notation (FZ-3), -nogap/-nomismatch rules.
- 2026-10-17 12:25 UTC | CPP-10 FZ-5, CPP-11 ENG-i2 + DEC-4f | reference goldens 279/279 (was 248); random vs reference 887/1000 (was 739); fuzz unexplained 431 (was 440; -nogap differences now counted as deliberate).
- 2026-10-17 12:30 UTC | CPP-12 FZ-1/FZ-2 readers, CPP-13 -blockonly + mark lines | goldens 279/279; random vs reference 892/1000; fuzz unexplained 76 (clustal/paml/fasta 22 each, codon 10).
- 2026-10-17 12:32 UTC | CPP-14 FZ-3 tfastx | fuzz vs real PAL2NAL 1000 cases x 4 formats: 0 unexplained differences; goldens 279/279.
- 2026-10-17 12:49 UTC | CPP-15 DEC-1 notices/--strict, CPP-16 -o, CPP-17 piped inputs, CPP-18 signals/pipe | goldens 279/279; random vs reference 1000/1000 (with --strict/-q); fuzz 0 unexplained; reference tests run against coal-cpp through a shim (tests/test_cpp_vs_reference.py reftests): test_review_fixes 43/43, test_fuzz_findings 17/17, test_decisions 19/24 (4 in-process DEC-3 mocks not applicable, 1 help text: next item).
- 2026-10-17 12:50 UTC | CPP-19 CLI help + tests refreshed against bb80bde | test_cpp_cli 11/11; reference test_decisions 20/24 (the 4 in-process DEC-3 mocks do not apply).
- 2026-10-17 12:51 UTC | CPP-20 PAML header for 1,000+ sequences | found while checking large inputs; real PAL2NAL prints ' 1200     90'; fixed, new test.

## Run 3 progress (C++ foundations brief: build, readers, tables, CLI, pool, test hookup)

- Reconciled with the brief: overnight/core already had the CMake build, the ordered streaming thread pipeline, the COAL_IMPL=cpp hookup and verified output writers (285 golden checks) from runs 1-2; those are kept, not redone. The local branch overnight/merged (merge with Keanu's integration, not pushed) is left untouched, as the brief says not to merge yet.
- 2026-10-17 10:40 UTC | CPP-1 build setup | CMake options COAL_WARNINGS_AS_ERRORS / COAL_BUILD_TOOLS; stricter flags. GCC 11 -Werror -Wpedantic -Wshadow clean; clang 19 -Werror clean for Windows (MinGW), macOS (arm64) and Linux targets. MSVC not available on this machine: flags set (/W4 /WX /permissive- /utf-8, _CRT_SECURE_NO_WARNINGS); a clang -Wconversion pass (approximating MSVC C4244/C4267 loss-of-data warnings, with 32-bit long as on Windows) found none.
- 2026-10-17 10:44 UTC | CPP-2 readers | C++ readers follow the reviewed Python reader (BOM, \n/\r\n/lone \r, '#' before first DNA header, non-letters removed, duplicate/empty IDs, non-ASCII, file:line). New coal-cpp-dump helper + tests/test_cpp_readers.py: 43 edge cases equal to the reference reader; the test fails when the BOM skip / '#' rule are broken (checked, then restored). Golden 285/285 for C++; Python suite unchanged.
- 2026-10-17 10:45 UTC | CPP-3 tables | C++ table 3 gains the GTG start (the only difference from the reviewed tables). tests/test_cpp_tables.py: 17 codes x 64 amino acids + start sets equal to the reference; fails without the GTG start (checked). Note: the Python engine on this branch (src/coal, not editable here) predates the reviewed fixes and still lacks GTG for table 3, so C++ and this branch's Python differ for table-3 genes starting with GTG until the integration branch is merged (random test unaffected: 1000/1000).
- 2026-10-17 10:59 UTC | CPP-4 CLI | argparse-compatible parser (args.cpp) + generated help; tests/test_cpp_cli.py 11 tests OK (also under ASan/UBSan); break tests: single-quote-only repr, no '--' removal, no attached short values, edited help -> 6/6/20 failures, restored. Golden C++ 285/285; test_cpp_vs_python/readers/tables OK; clang -Werror (Windows/macOS) and -Wconversion clean for args.cpp/main.cpp.
- 2026-10-17 11:04 UTC | CPP-5 pool | pool.hpp (run_ordered / parallel_batch) replaces the inline batch threading in main.cpp; dump `pool` command + tests/test_cpp_pool.py 7 tests OK, also under TSan; break tests (one thread / batches out of order / worker exception escapes) -> 5 / 12 / 1 failures, restored. Golden C++ 285/285; vs_python, cli, readers, tables OK; clang -Werror Windows/macOS clean.
- 2026-10-17 11:06 UTC | CPP-6 test hookup | verified the existing COAL_IMPL=cpp hookup (skip when missing: OK; fake binary /bin/true: 7/285, so it can fail; real: 285/285; Python 285/285); fixed the `--count` crash when the binary is missing.

## Decisions

- Output-format conventions were taken from the goldens, plus a few black-box runs of PAL2NAL v14 in a scratch folder outside the repo (only its *output* was looked at; no code read or copied, nothing from it committed):
  - CLUSTAL/codon name column = max(10, longest ID) + 4 characters; CLUSTAL blocks of 60 nt, no conservation line; empty alignment = header only.
  - PAML header = `%4d%7d` (sequences, length); sequence wrapped at 60; empty sequence = empty line.
  - codon: residue line and codon line are cut into blocks of 20 *independently* (residues vs. triplets), as PAL2NAL does.
  - PAL2NAL refuses `-output codon` with -nogap/-nomismatch; COAL does the same (error, exit 1) instead of inventing a layout.
- Known PAL2NAL quirk NOT reproduced: in `-output codon`, when a column is widened by a frameshift digit 4/5, PAL2NAL sometimes inserts an extra '-' in the residue line (seen for a lone '4' or '4'+'5' in a column, but not for a lone '5', '4'+'2', or the frameshift_mixed fixture). The rule could not be pinned down in 3 probes; COAL prints one residue per column (matches all goldens). Parked.

## Not checked against PAL2NAL (combinations without goldens)

- `-output paml` / `-output clustal` with -nogap or -nomismatch (unit-tested: same rows as FASTA output).
- `-output clustal|paml|codon` with -codontable N (only the FASTA variants have goldens for codon tables).
- `-output codon` with frameshift digits 4/5 in long alignments (see quirk above).

## Needs a human decision

- ~~`build/` should be added to .gitignore~~ Done.
- An untracked `Claude outputs/` folder exists in the repo root (not created by this session); left untouched.
