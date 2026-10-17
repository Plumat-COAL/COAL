# Overnight log: core (branch overnight/core)

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

## Progress log

- 2026-10-17 01:17 UTC | task 0 | preflight OK: cmake 4.4.4 + g++ 11.4, C++17 + std::thread hello-world built and ran (2 hardware threads). Branched overnight/core from main 31cf3d0 (pull: already up to date).
- 2026-10-17 01:17 UTC | task 1 | Python 279/279 golden checks (all variants: fasta, clustal, paml, codon, nogap, nomismatch, codon tables), C++ n/a | added -output clustal/paml/codon; test_differential checks every variant; -j1 == -j2/-j3 tests for every format.
- 2026-10-17 01:25 UTC | task 2 | Python 279/279, C++ 279/279 golden checks | C++17 core in cpp/src (tables, core, seqio, main) + root CMakeLists.txt (coal_core static lib + coal-cpp exe). Whole port compiled warning-free and passed every golden check on the first build, so there were no intermediate counts to commit. test_differential runs either engine (COAL_IMPL=python|cpp, COAL_CPP_BIN) and has a --count mode.
- 2026-10-17 01:27 UTC | task 2 | Python 279/279, C++ 279/279 | tests/test_cpp_vs_python.py: 1,000 random inputs (seed 20261017; Python outcomes 301 ok / 647 with warnings / 52 errors) -> 1000/1000 identical stdout + full stderr + exit code, ~7 s; plus -j1/-j2/-j3 identical on a 3.4M-nt, 2-batch dataset for all formats, warnings in input order.
- 2026-10-17 01:31 UTC | task 2 | Python 279/279, C++ 279/279, random 1000/1000 | C++ speed-ups: table lookups hoisted out of hot loops, sequential (not seek-per-row) reads of staged rows, memchr line scanning in the DNA index. 10k x 500 on the VM's local disk: C++ 0.09 s (2 threads) / 0.12 s (1 thread) vs Python 0.59 s. run_bench.py: --impl python|cpp|both, --cpp-bin.
- 2026-10-17 01:34 UTC | task 2 done | Python 279/279, C++ 279/279, random 1000/1000 | expand() fast path, token reserve. Measured (2-core VM, local disk, medians): 10k x 500 Python 0.593 s / C++ 0.097 s (6.1x); 100 x 10k 0.149 / 0.023 (6.4x); 100k x 500 3.559 / 0.881 (4.0x); C++ peak < 30 MB. Results in benchmarks/results/cpp_vs_python_linux.md (+ csv/txt).
- 2026-10-17 01:40 UTC | tasks 3+4 | Python 279/279, C++ 279/279, random 1000/1000, all suites OK | src/coal/report.py: --report FILE.json / --json (needs -o), schema coal.report/1 exactly as specified plus extra fields (sequences[].index, summary.not_processed, summary_text, error{code,text,detail}, hints[].checked/fit/fit_current); report written on failed runs too (plain-English reason per error code). End-of-run 'coal: summary:' lines + genetic-code hint (>20% mismatching/failed -> re-test up to 200 of them with all 17 codes; suggest only if it fits >=80% and >=20 points more than the current code; ties -> lowest table number). The summary is collected on every run with counters only (per-sequence entries only kept when a report is requested). -q suppresses summary lines. tests/test_report.py (10 tests), tests/report_example.json (real run: 6 mito genes with code 1 -> hint code 2). C++ does not print summary lines yet; test_cpp_vs_python ignores 'coal: summary:' lines (allowed by the task).

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

- `build/` should be added to .gitignore (I may not edit .gitignore).
- An untracked `Claude outputs/` folder exists in the repo root (not created by this session); left untouched.
