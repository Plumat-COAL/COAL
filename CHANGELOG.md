# Changelog

## Unreleased

- **GUI: the alignment view's residue-class colours are readable for colour-blind users.** The seven classes keep
  their Okabe-Ito hue families but now sit on an even lightness ladder, so they stay apart under protanopia,
  deuteranopia and tritanopia and in greyscale (the pastel tints of equal lightness merged there; verification
  2026-10-24, V-17); the text stays at 5:1 or better on every tint, in the light and the dark theme. A test holds
  the simulated distances. The letters, outlines and legend carried the meaning before and still do.

Found by the full verification of 2026-10-24 (`verify/REPORT.md` on the branch `verify/full-2026-10-24`: 16 minor,
usability, docs and test findings, no wrong output):

- **`pal2nal` without arguments is a usage error:** the usage text goes to standard error with exit 1 (as the
  drop-in's help says and as PAL2NAL's own usage goes to stderr), not to standard output with exit 0, so a pipeline
  step that forgot its files no longer gets the usage text as its "alignment". `pal2nal -h` is unchanged (stdout, 0).
- **GUI: a run that stops before pairing says why in the banner** (the engine's own sentence, e.g. "control character:
  sequence files must be plain text (is this a binary or word-processor file?)", file names only); it used to say
  only "COAL stopped: no alignment was written" with the reason inside the collapsed raw-messages box.
- **GUI: the genetic-code tip on a large input says it tested a sample** ("50 of the 50 genes COAL tested (a sample
  of the 100,000)") instead of "50 of 50 genes".
- **GUI accessibility:** the "Methods paragraph" box is announced by its label (its copy button, which sat inside the
  label and took the label for itself, is now a separate "Copy methods paragraph" button), and Gradio's hidden
  duplicate tab bar no longer holds focusable buttons (axe `aria-hidden-focus`).
- **Docs:** the desktop README's privacy line says none of *your data* leaves the computer and names the one helper
  script Gradio's page requests from a CDN (not needed; the page works with it blocked); the installer name is the
  version-free `COAL-Setup-windows-x64.exe`; `docs/usage.md` says a `kill -9` of a run writing to standard output can
  leave `.coal_<random>.tmp` in the temp folder; `tests/README.md` says which of the "17 local codon alignments" are
  in the repository (12) and which are local-only (5); SKILL.md's "after the engine update" and "not yet verified on
  Windows" sections are current (both engines verified on Windows 2026-10-24); `demo/mcp_demo.md` regenerated
  (`codon_align` keeps going past a broken gene and returns `completed: true` with `exit_code: 2`).
- **Tests:** `test_docs_explainer` finds Playwright's Chromium on every platform (it looked for a Linux-only path and
  was always skipped elsewhere); `test_cpp_vs_reference` copies instead of symlinking where symlinks need a privilege
  (a stock Windows account); `test_resume_with_workers` ignores a venv launcher's same-argv child on Windows; new
  tests for the GUI's download buttons (the one gap found by mutation testing: 12 of 15 planted bugs were caught),
  the stop reason, the sampled tip and the methods label. The PAL2NAL fuzz harness no longer reports a difference
  for `-codontable 10` genes that both programs place by position 1 or 2 nt off the generator's frame (the
  neutralised TGA then straddles two codons of the alignment grid; seeds 20261024 and 20261025 each had one such
  case, with COAL's output equal to PAL2NAL's on the original input). `test_cpp_vs_reference` is skipped on
  Windows with the reason: its shim runs coal-cpp as a child process there, so the reference tests that interrupt
  a run watch the wrong process (before, the test errored on the symlink privilege before getting that far).

Found by the stress test of 2026-10-23 (`docs/stress-test-2026-10-23.md`), both engines unless said:

- **Two runs writing the same `-o` file at the same time are refused** instead of sharing one `<out>.coal-partial`
  file, which could end a run with exit 0 and an empty or another run's output. A run with `-o` keeps
  `<out>.coal-lock` while it runs; a lock left by a killed run is replaced once its process is gone (a lock
  file that is still empty is one being created by another run, not a stale one: the second iteration of the
  stress test caught two runs reading each other's empty lock and both going ahead).
- **Thousands of DNA files:** at most about as many files as the open-file limit allows stay open (least recently
  used closed first, the limit raised towards its hard limit when needed) and their read buffers together take
  about 64 MB at most. 1,000 files failed with "too many open files" in the Python engine with workers, 60,000 open
  files took 8 GB (129 KB each on Python 3.14), and more than 65,535 files overflowed a 16-bit file number: a
  traceback in the Python engine, a silent wrap in coal-cpp.
- **coal-cpp: the exact-fit search is linear in the DNA** (a string search over the translated frames) instead of
  trying every start codon by codon: a 200,000-residue low-complexity gene behind a near-copy of itself took 138 s
  (an hour at 1,000,000 residues), now well under a second, as in the Python engine. Same answer everywhere
  (the differential tests and the engine fuzz).
- **A lone `\r` ends a line anywhere in a DNA file** (universal newlines, as the protein reader): the style used to
  be decided from the first 64 kB, so a file whose line ends change part-way lost every record after the first lone
  `\r` beyond that and failed with a misleading count error. PAL2NAL v14 reads such files.
- **`nohup coal ...` survives the terminal closing:** an already ignored SIGHUP stays ignored (the handler used to
  replace nohup's setting; the run died with exit 129).
- **Out of memory in the main process is one `coal: error: out of memory ...` line** (exit 1, report code
  `out_of_memory`) instead of a traceback (a worker dying of it was already handled).
- **App: `python -m coal_app` no longer dies on a non-UTF-8 console** (Gradio's MCP banner emoji on a cp1252
  Windows console or pipe).
- **`-q` without `--report`/`--json` keeps no per-sequence entries** and writes none to the `--resume` entries
  file (a million genes of unrelated DNA wrote 55 GB of it), so memory and disk stay flat on data where every
  codon mismatches (the summary lines, the hint and the report, which read them, are not produced under `-q`).
- **Memory no longer grows with the number of per-codon warnings, in either engine** (finding 9 of the stress
  test, both halves). The run report keeps, per sequence, one byte for its kind, one for a few flags and the
  codon of its first definite mismatch (plus the first flagged and the first excluded gene whole, for the
  summary lines, and where up to 200 failed genes break), instead of every per-sequence entry with every message
  until the end of the run; the entries go to the `--resume` entries file next to the output as the run goes
  (or to a temporary file when a report is wanted without `-o`; nothing with `-q` and no report), and `--report`
  / `--json` read them back in order. A resumed run rebuilds its state from that file one line at a time (the
  checkpoint check no longer loads the whole file either). The Python engine's worker processes hand back each
  sequence's mismatch columns and codons instead of its warning strings (about 2 KB against 35 KB for a
  400-codon gene of unrelated DNA, per gene of the batch in flight in every worker and in the main process), and
  the main process formats a warning only when it prints it; the strings, the report, stderr, the exit status
  and `--resume` are unchanged byte for byte (`tests/golden_memory`: recorded from the previous `main`, both
  engines). Measured on 400-codon genes of unrelated DNA, every codon a warning (peak memory as the summed PSS
  over the process tree; 16 threads): 20,000 genes in default mode, coal-cpp 1,172 MB before and 103 MB after,
  the Python engine 2,822 MB and 45 MB (`-j 1`), 3,226 MB and 338 MB (`-j 16`, 7 workers); 100,000 genes,
  coal-cpp 5,654 MB and 118 MB, the Python engine 13,761 MB and 69 MB (`-j 1`), 14,647 MB and 885 MB (`-j 16`,
  16 workers: 545 MB of that is what the same input in clean form costs with 17 processes, the rest the batches
  in flight). The default mode is also faster on such
  data (100,000 genes, `-j 1`: 213 s to 158 s in Python, 185 s to 161 s in coal-cpp) since no entry store grows;
  the one cost is the main process formatting the strings the workers formatted in parallel before, which only
  shows when stderr takes most of a run. A run with `-o` on such data writes the entries file next to the
  output while it runs (about 55 KB per 400-codon gene of unrelated DNA, removed at the end).
- **The genetic-code hint no longer re-tests genes of more than 100,000 residues** (28 s of a 30 s run on one
  1,000,000-residue gene with one wrong codon; `-q` already skipped it).
- **Messages name at most three input files** ("a.fa, b.fa, c.fa and 997 more files") instead of every path.
- `tests/tools/stress_*.py`: the stress-test generator, harness and helpers; `tests/test_stress_findings.py`.

## 0.5.0

- **Run report: a new top-level `engine` field** (`coal.report/1`, next to `source`): `"coal-python"` or
  `"coal-cpp"`, the same strings as the `--resume` checkpoint's `engine`, so a pipeline or the app can tell which
  engine wrote a report (also when `coal` handed a large input to `coal-cpp`). The per-job reports of a
  `coal --batch` report (`coal.batch/1`) carry it too. Apart from `engine` and `runtime_seconds` the two engines'
  reports stay identical (#18).
- **Changed: `summary.not_checked` of a run that stopped at an error now counts only the sequences the run never
  reached**, the number the `coal: summary:` line "N sequences were not checked" gives; the sequences checked
  before the error are counted as what COAL found (`ok`, `warning`, `failed`), so `ok + warning + failed +
  not_checked == total`. It used to count every sequence but the failed one (e.g. `examples/one_broken_gene`
  without `--keep-going`: `ok 0, failed 1, not_checked 19`; now `ok 13, failed 1, not_checked 6`). The per-sequence
  entries are unchanged (`not_checked` with `not_written` before the error, `not_reached` after it). Both engines;
  documented in `docs/report-schema.md` (#18).
- **`pip install coal-align` brings the C++ engine** (issue #16): each release now has one wheel per platform with
  `coal-cpp` inside, at `coal/bin/coal-cpp` in the installed package: Linux x86_64 and aarch64 (`manylinux_2_28`,
  glibc 2.28+, `libstdc++` linked statically), macOS 11+ (one `universal2` wheel for Apple silicon and Intel) and
  Windows x64 (static CRT). pip picks it by itself from PyPI; from a GitHub Release, install the wheel file for your
  platform. The pure `py3-none-any` wheel and the sdist stay as the fallback everywhere else. `coal` finds the bundled
  binary with no settings: the lookup order is now `COAL_CPP`, then `coal/bin` inside the package, then next to
  Python, then `PATH` (docs/cpp-engine.md "Where coal looks for coal-cpp"). New `coal --cpp-path` prints the
  `coal-cpp` that `coal` would use (exit 1 with the reason when there is none or it is another version), for
  pipelines that call the binary directly. The wheels are made by `packaging/wheels/build_platform_wheel.py` (no
  cibuildwheel, no compiled extension) in `release.yml`'s new `cpp-wheels` job and checked there (auditwheel, lipo, a
  fresh-venv install whose output equals the Python engine's); `tests/test_wheels.py` does the same with a local build.
- **Nightly fuzz workflow** (`.github/workflows/nightly.yml`, 03:17 UTC, or by hand with a seed): the fuzz against
  the real PAL2NAL v14 with a fresh seed each night, once per engine (`COAL_FUZZ_ENGINE=cpp` in
  `tests/fuzz_vs_pal2nal.py` runs `coal-cpp`), the new `tests/tools/engine_fuzz.py` (coal-cpp vs the Python engine
  on malformed inputs and every option; Linux and macOS) and the C++ tests plus the generator under
  AddressSanitizer and UndefinedBehaviorSanitizer (`COAL_SKIP_MEMORY_TESTS=1` skips the peak-memory test there).
  Differing cases are uploaded as artifacts; CONTRIBUTING.md "Nightly checks" says how to replay them.
- **coal-cpp: no read past the end of a short gene's DNA** when the expected position of an anchor (after a
  frameshift digit or a mismatch) lies beyond the sequence end; found by the first sanitizer run. Output and exit
  codes are unchanged (the out-of-bounds read could never have produced a match).
- **MCP server on mcp 1.x (the Windows app's) is quiet and names its version:** `python -m coal_app --mcp` /
  `COAL.exe -m coal_app --mcp` with the `mcp` 1.x that Gradio's extra installs reported the mcp library's version
  ("1.10.1") as `serverInfo.version` and logged about 37 KB of INFO lines to standard error per session (into
  Claude Desktop's MCP log). It now reports COAL's version, as on mcp 2.x, and logs warnings and errors only. The
  message for a missing `mcp` package no longer claims that version 2 is required (1.x works too).
- **Packaging:** the source distributions of `coal-align` and `coal-align-app` no longer ship `tests/` (the tests
  need the repository: helpers, fixtures, examples, coal-cpp), and the `coal-align-app` wheel and sdist now carry
  the MIT `LICENSE` (the app had none).
- **App: a settings file that is not text no longer stops COAL from starting.** `settings.json` in another encoding
  (or binary garbage) gave a traceback at start-up; it now counts as a damaged file: the defaults are used and the
  Settings tab says so, as for a half-written file.
- **App: an unusable results folder is one plain sentence, not a traceback.** With `COAL_RUNS_DIR` pointing where no
  folder can be made (a path under a missing directory, an existing file, no permission), the GUI shows "Could not
  create the results folder <folder>: <reason>" in the banner, the MCP tool `codon_align` raises a tool error with
  that sentence, and `coal-gui` / `python -m coal_app` print it and exit 1 instead of failing in `os.makedirs`.

- **Device inputs are refused** (robustness review): `coal /dev/zero nuc.fa` or `coal - nuc.fa </dev/zero` used to
  copy the device to a temporary file without end, filling the disk silently. Both engines now stop before reading
  with `coal: error: cannot read /dev/zero: it is a device, not a file` (exit 1; `standard input` for `-`). A
  terminal on standard input and `/dev/null` behave as before; pipes, FIFOs and `<(...)` are unchanged.
- **`--report PATH` no longer replaces what is at PATH** (robustness review): it used to rename a temp file over
  anything there, so `--report nuc.fa` destroyed an input, `--report out.fa` the output, and a symlink, FIFO or
  device node (`/dev/null` became a regular file!) was replaced. Both engines now give the report path the rules of
  `-o`: an input or the output file is a usage error before anything runs (`coal: error: the report path r.json is
  an input file` / `is the output file`, exit 2; also for a job line's `--report` and the batch report), symlinks
  are written through, and FIFOs, devices and hard-linked files are written in place. A report into a folder that
  does not exist is a usage error too (`the report path results/r.json: folder /data/results does not exist`, exit
  2) instead of a full run that ends by losing its report, and the report's temp file is removed on SIGTERM /
  SIGHUP / Ctrl-C while it is written, in both engines.
- **A killed worker no longer hangs the Python engine** (robustness review): when the out-of-memory killer (or
  `kill -9`) took one worker of a `-j` run, the main process could block forever in the pool's shutdown. It now
  ends the pool without waiting and stops with `coal: error: a worker process died (out of memory?)` (exit 1, report
  code `worker_died`, no partial output, no temp files). Workers also let SIGTERM end them (the inherited handler
  used to be swallowed inside the work) and **die with the main process** (`PR_SET_PDEATHSIG` on Linux, a parent
  watchdog everywhere), so a `kill -9` of `coal` no longer leaves workers holding its stdout open and a pipeline
  reading it waiting for EOF forever.
- **Python 3.14 on Linux: every multi-worker run failed** with `a worker process died (out of memory?)` after
  half a second. 3.14 made `forkserver` the default start method there, so the workers' parent is the fork server,
  not `coal`, and the watchdog above took that for a dead parent. The workers now wait for EOF on a pipe whose only
  write end `coal` holds, whatever the start method; that also fixes **Windows, where a killed `coal` left its worker
  processes running** (there `getppid()` keeps returning the dead parent's pid). `COAL_MIN_NT_PER_JOB` (default 4,000,000)
  can now be set in the environment, as `COAL_BATCH_NT`, so tests can force worker processes on a small input; CI
  runs the suite on 3.14 (Linux).
- **Worker processes end as soon as the main process fails** (Python engine): after an error on the writing side
  of a `-j` run (a full disk, a file-size limit) the pool was only ended when the cyclic garbage collector got to it
  (the suspended result generator and an exception held in a local kept it alive), so terminated workers and pool
  threads could outlive the error by seconds. The generator is now closed on the way out and the pool ended at once.
- **Checkpoints are checked before `--resume` trusts them** (robustness review): an edited or corrupted
  `<out>.coal-checkpoint.json` (`completed` beyond the input or negative, a string or float where an integer belongs,
  a list at the top level, a missing key, a deleted entries file, a partial file shorter than recorded, a truncated
  or deeply nested file, a folder at the checkpoint's path) made both engines write a truncated output with exit 0,
  or crash with a traceback / raw C++ text. Both engines now check every field's type and range and the partial
  files' sizes, treat any problem as the existing "unreadable checkpoint" notice and start a normal run; a
  checkpoint that cannot be written (`cannot write the checkpoint ...: is a directory`) is a notice and the run
  goes on without one; and `coal-cpp`'s JSON parser has a nesting limit (64) instead of overflowing the stack.
  The `drop` and `bad` column lists are bounded by the run's protein and codon columns too: a crafted
  `"bad": [4000000000000]` made `coal-cpp` fail with a raw `std::bad_alloc`; it is an unreadable checkpoint now, in
  both engines.
- **`--resume` of a `--keep-going` run in fasta or paml format** (Python engine): a checkpoint written after a
  skipped sequence was rejected as "unreadable" (the skipped sequence's index was recorded for clustal / codon
  only, and the checkpoint check requires one index per skipped sequence), so the run started over. The index is
  now recorded in every format, as `coal-cpp` does, and the resumed output is byte-identical to an uninterrupted run.
- **A checkpoint written by the other engine** (robustness review): `coal --resume` now picks the engine that wrote
  the checkpoint (so a `coal-cpp` installed or removed between the run and the resume does not get in the way),
  and when an engine is asked to resume the other one's run it says so instead of blaming the inputs: `out.fa was
  interrupted by the Python engine (coal-python), which must resume it: run again with COAL_ENGINE=python, or start
  a new run without --resume`.
- **`-o /dev/fd/1` and `-o /dev/stdout` work when stdout is a pipe** (robustness review): both engines resolved
  the path to `pipe:[N]` and failed with `cannot write /dev/fd/1: No such file or directory`. An existing target
  that is not a regular file is now written through before any path resolution.
- **`--report /dev/fd/1` and `--report /dev/stdout` work when stdout is a pipe**, in both engines, for the same
  reason (the report failed with `could not write the report /dev/fd/1: No such file or directory` after the
  whole run).
- **`--json | head` is a silent exit 141** (robustness review) in both engines (the Python engine printed a
  traceback and exited 1, `coal-cpp` exited 0), also for `--batch ... --json`. `coal-cpp` now streams the `--json`
  report to standard output one sequence at a time instead of building it whole in memory first (the bytes are
  the same).
- **`coal-cpp` falls back when `TMPDIR` is unusable** (robustness review): a `TMPDIR` naming a missing folder or a
  file made it fail with `filesystem error: temp_directory_path`; it now searches `TMPDIR`, `TEMP`, `TMP`, `/tmp`,
  `/var/tmp`, `/usr/tmp` and the current folder as Python's `tempfile` does, so both engines behave the same.
- **Disk-full errors name the file** (robustness review): the Python engine printed `[Errno 28] No space left on
  device` (and could hang at exit after terminating its workers), `coal-cpp` `could not write the output (disk
  full?)` (and, with the rows going through a full temp folder, exit 0 with an empty output). Both now say
  `cannot write out.fa: no space left on device`, `cannot write a temporary file in /tmp: ...` or `cannot copy
  standard input to a temporary file in /tmp: ...`, exit 1, and leave no temp file behind.
- **Python API checks its arguments** (robustness review): `coal.codon_align` raises `CoalError` with
  `code="duplicate_id"` for a repeated protein or DNA ID (it used to pair them silently; the command line refuses
  the file), and both `codon_align` and `coal.align_files` raise `TypeError` with a message naming the argument for
  a dict, generator, bare string, `None`, a non-`str` ID or a pair of the wrong length, instead of a `KeyError` or
  `AttributeError` from inside the engine. `align_files(pep, "nuc.fa")` with one DNA path as a string or
  `os.PathLike` now means that one file (the string used to be iterated character by character). `nuc_paths` is
  one path or a `list` / `tuple` of paths: a generator (or any other iterator) of paths, which used to work, is a
  `TypeError` now; pass `list(...)`.
- **Smaller robustness fixes** (both engines, identical messages): a NUL or other control byte (except tab, LF,
  VT, FF, CR) anywhere in a protein or DNA file is `<file>:<line>: control character: sequence files must be plain
  text (is this a binary or word-processor file?)` (exit 1, code `input_format`) instead of a confusing row-length
  error, a mangled ID or silently dropped bytes; IDs longer than 200 characters are shortened on stderr to their
  first 100 and last 20 characters (the output and the report keep them whole); `-o -` means standard output; the
  Python engine reports a DNA file that shrank between indexing and reading as `could not read the DNA of record N
  (the file changed after it was read, or cannot be read)`, as `coal-cpp` does, instead of reading another
  record's bytes; an empty DNA file among several is documented as accepted.
- **Mismatching long genes no longer take quadratic time.** When a protein does not fit its DNA exactly, both engines
  locate PAL2NAL's 10-residue chunks in the DNA; every chunk was a scan from the start of the DNA (coal-cpp; the
  Python engine scanned for every chunk that was not where the previous ones put it), so a 35,000-codon gene with one
  wrong codon cost about 10 s per gene with coal-cpp (and the genetic-code hint repeated that under 17 other codes).
  The DNA is now translated once per frame and its 10-residue windows indexed, so a chunk of ordinary residues is a
  lookup; chunks with a frameshift digit, an X or the first Met keep the old scan. Same chunks found, same output and
  messages: 10 genes × 35,000 codons with one wrong codon each, 104 s → 3.0 s (coal-cpp), 27 s → 16 s (Python);
  exactly fitting genes are unchanged. The Python engine's fallback is also about 1.5× faster per short gene.
- **coal-cpp uses far less memory on runs with many warnings or with `--report`.** It kept a JSON tree for every
  sequence with a warning (about 4 KB each: 1 GB at 300,000 mismatching sequences) and, for `--report`/`--json`, built
  every sequence's entry before writing (3 KB each: 3.2 GB at 1,000,000 sequences). It now keeps each non-ok entry as
  its compact JSON line (about 0.6 KB; the Python engine's dict entries take 1.4 KB) and streams the report one
  sequence at a time:
  300,000 mismatching sequences 994 → 266 MB; 1,000,000 clean sequences with `--report` 3.2 GB → 204 MB. The
  Python engine streams its report the same way and no longer builds a view of every sequence for the genetic-code
  hint (1,000,000 sequences with `--report`, one process: 635 → 367 MB). The report JSON is unchanged.
- **Docs:** `docs/usage.md` gives the measured memory per sequence of both engines (about 350 bytes with the Python
  engine after the change above, 190 with coal-cpp; the old figures, 300 bytes and 122/139/145 MB, matched coal-cpp
  only), where each engine puts its temporary files and how much disk the clustal/codon formats need, that
  `-q --progress` still prints progress, and that the checkpoint cadence is best-effort under I/O contention.

- **Less memory with many workers:** COAL's own processes now limit numpy's BLAS thread pool to one thread
  (`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, unless already set). COAL uses no BLAS, but every
  worker that imported numpy started a pool sized to the machine. On a 20-thread Windows machine each process that
  imported numpy committed about 630 MB (9-11 GB at the default `-j`); the BLAS pool is the likely cause, to be
  confirmed by re-measuring. Programs that use COAL as a library are not affected.
- **Faster start-up** for runs on small files (one COAL call per gene family in a pipeline): about 40 ms less per call
  (104 -> 64 ms on Linux), by not importing `dataclasses`, `subprocess`, `hashlib`, `json` and `shlex` unless needed.

- **`-` reads standard input** (`zcat nuc.fa.gz | coal pep.aln -`), for the protein alignment or one DNA file, in
  `coal`, `coal-cpp` and the `pal2nal` drop-in, also on Windows (no `/dev/stdin` needed). The data is read as bytes
  into a temporary file, as for other pipes, so the output is byte-identical to giving the file's path; messages name
  it `standard input`; a second `-` is a usage error. The run report records `"path": "-"` with the sha256 of the data
  read (the report of a run reading `/dev/stdin` or another pipe now hashes the data too, in both engines, instead of
  the Python engine's hash of an empty re-read).
- **C++ library for embedding** (e.g. in DIAMOND): `cmake --install` puts `coal-cpp`, `libcoal_core.a`, the headers
  (`<coal/core.hpp>`, `<coal/tables.hpp>`, `<coal/seqio.hpp>`, `<coal/version.hpp>`) and a CMake package in place;
  `find_package(coal CONFIG)` or `add_subdirectory` / FetchContent give `coal::core`. As a subproject only the library
  is built (`COAL_BUILD_CLI`, `COAL_BUILD_TOOLS`, `COAL_INSTALL` to add more). Example: `examples/cpp-embed/`; docs:
  docs/cpp-engine.md "Embedding the library".
- **orthologr stand-in takes every PAL2NAL option:** `packaging/orthologr/pal2nal.pl` sent every command line to
  `coal-cpp` when installed, which rejects PAL2NAL-only options such as `-nostderr` and prints `coal`'s help for
  `-h`. It now runs `coal-cpp` only for the options it understands (`-output`, `-codontable`, `-nogap`,
  `-nomismatch`, `-blockonly` and the files) and COAL's `pal2nal` drop-in for anything else; `COAL_PAL2NAL` still
  overrides.
- **orthologr:** `packaging/orthologr/pal2nal.pl` is a stand-in for the `pal2nal.pl` bundled in the orthologr R
  package; copied over it, orthologr's `codon_aln()` / `dNdS()` run COAL (coal-cpp when installed) with their
  unchanged command line and get the same output. Notes on what to expect and on a proper `tool = "coal"` option
  for orthologr: `packaging/orthologr/README.md`.

- **Python API: every error is a `CoalError`** (issue #18). `coal.codon_align` and `coal.align_files` with an
  unsupported `table` now raise `CoalError` with `code="unknown_table"` (the message lists the supported NCBI codes)
  instead of a bare `ValueError`.
- **Python API: `coal.seqio.InputFormatError` is a `CoalError`** with `code="input_format"` (issue #18), so one
  `except CoalError` catches a file that is not FASTA or CLUSTAL too; it stays a `ValueError` for code that catches
  that. The `coal` command's messages, exit code and report are unchanged.
- **`pal2nal | head` no longer prints a traceback** (issue #18): the drop-in's usage text into a closed pipe is a
  silent exit 141, as for `coal` and for the alignment output.
- **Batch mode (`coal --batch JOBS`, #15):** thousands of small conversions in one process, one `coal` command line
  per line of JOBS (`-` = stdin), the command line's options as the defaults of every job, `--batch-out-dir DIR` for
  jobs without `-o`, a `coal.batch/1` report with `--report`/`--json`, per-job stderr prefixes
  (`coal: [k/n PEP] warning: ...`) and a `coal: batch:` summary line. A failing job does not stop the batch (exit 1
  at the end; 2 if jobs only left sequences out; mistakes in JOBS are reported with their line number before
  anything runs). Every output is byte-identical to the single run; about 40× less wall time on tiny jobs. Python
  only: `coal-cpp --batch` is a clear error and batches are never dispatched to it. `--resume` and cross-job
  parallelism are not supported in a batch yet.
- **The app no longer waits forever for a stopped engine:** when an engine run hits the app's time limit (or the
  demo pipeline's), it is sent SIGTERM and, after the grace period, SIGKILL; the kill now takes the engine's whole
  process group, so its worker processes die with it. Before, a run that outlived the grace period on a slow machine
  was killed alone, its orphaned workers kept its output pipes open and the app waited for them indefinitely (the
  CI app suite hung until the job timeout). Windows already killed the tree. Test commands that start the engine
  now give it no stdin and always carry a timeout.
- **`coal-cpp` mirrors Python 3.13's argparse** (it mirrored 3.10's): `-qx` is "unrecognized arguments: -x" rather
  than an ignored explicit argument, `-hx` prints help, `-cod=2` is accepted, invalid choices are quoted
  (`invalid choice: '7' (choose from '1', '2', ...)`), and only the `--` marker is dropped from the positionals (a
  second `--` is a file name). `tests/test_cpp_cli.py` compares the two command lines live on an interpreter
  whose argparse behaves that way (probed, not keyed on the version) and skips elsewhere with the differences named.
- **Documented for 1.0:** `docs/python-api.md` (the public Python API: `coal.codon_align`, the new re-export
  `coal.align_files` = the `coal` command as a function, `coal.CoalError` with its `code`/`seq_id`/`index`, and the
  semantic-versioning promise for command-line flags, output bytes, exit codes, the report schema and these
  functions), `docs/report-schema.md` with a JSON Schema (`docs/report.schema.json`, draft 2020-12) for the run report
  `coal.report/1` (every field, status, pairing mode, message and hint code; additive changes only within `/1`),
  `CONTRIBUTING.md`, GitHub issue forms and a pull-request template. New tests run the documented examples
  (`tests/test_python_api.py`) and validate freshly generated reports from both engines against the schema
  (`tests/test_report_schema.py`); `tests/report_example.json` is now a real `examples/one_broken_gene --keep-going`
  report. No change to COAL's behaviour or output.
- **orthologr:** `packaging/orthologr/pal2nal.pl` is a stand-in for the `pal2nal.pl` bundled in the orthologr R
  package; copied over it, orthologr's `codon_aln()` / `dNdS()` run COAL (coal-cpp when installed) with their
  unchanged command line and get the same output. Notes on what to expect and on a proper `tool = "coal"` option
  for orthologr: `packaging/orthologr/README.md`.
- **Repository tidied for outside users:** a short README with a section on using COAL in pipelines; the detailed
  documentation is in `docs/`. Hackathon logs and pitch material are removed (kept on the `hackathon-archive` tag),
  review-round markers are gone from comments, and engine test modules are named after what they test. No change to
  COAL's behaviour or output.
- **`coal-cpp` reports the same error as `coal` when a protein file has two defects.** A FASTA record is now
  handed back before the `>` line that ends it is checked, so a short row is reported before a later duplicate or
  empty header; and a CLUSTAL file's rows are parsed before its `#` mark line is read, so a malformed CLUSTAL line is
  reported before a later non-ASCII character. Found by differential fuzzing; `tests/test_cpp_error_precedence.py`.
- **Package renamed to `coal-align`** (`coal` on PyPI is another project). The `coal` command, the `coal` Python
  module and `pal2nal` are unchanged; release files are now `coal_align-<version>-py3-none-any.whl` and
  `coal_align-<version>.tar.gz` (`coal_src.tar.gz` stays).
- **The app is its own package, `coal-align-app`:** the web GUI, the MCP server, the viewer, the Windows app and the
  Hugging Face Space moved to `apps/coal-app/` (package `coal_app` unchanged; `coal-gui` and `python -m coal_app` work
  as before). `coal-align` is now the engine alone (`coal`, `pal2nal`, the examples) with no dependencies and no app
  code, so pipelines (OrthologR, DIAMOND workflows) install only that. The `gui` and `mcp` extras are on the app:
  `pip install "coal-align @ .../coal_src.tar.gz" "coal-align-app[gui] @ .../coal_app_src.tar.gz"`, or from a checkout
  `pip install -e . -e "apps/coal-app[gui]"`. Both packages always have the same version, and the app depends on that
  exact `coal-align`. Releases add `coal_align_app-<version>-py3-none-any.whl`, `coal_align_app-<version>.tar.gz`
  and the version-free `coal_app_src.tar.gz`. The Windows launchers are now `apps\coal-app\desktop\run_from_source.bat`
  and `apps\coal-app\desktop\setup_mcp.bat`. The app's tests (`apps/coal-app/tests/`) are grouped by feature; the
  checks of the docs, the website and CI moved to `tests/` (`test_docs_*.py`, `test_repo_ci.py`).

- **Bioconda recipe prepared** (`packaging/bioconda/`): `meta.yaml`, `build.sh` (coal-cpp via CMake into `bin/`, then the
  Python package), a fixture pair and a test that `coal` finds the installed `coal-cpp` and hands a 1 MB input to it.
  Not submitted yet; `packaging/bioconda/README.md` explains the sha256 step and the pull request to `bioconda-recipes`.

## 0.4.2

- **COAL in Claude Desktop, from the installer:** the Windows installer's last page has **Connect COAL to Claude
  Desktop** (ticked). Claude Desktop then starts COAL itself (`COAL.exe -m coal_app --mcp`) and can run it on your
  files; uninstalling removes the entry. `codon_align` can save the alignment where you ask (`save_as`), e.g. next to
  the input files.
- **C++ engine included:** the Windows app ships `coal-cpp.exe`, and `coal` hands inputs of 1 MB or more to it when it
  is found and is the same version (same output, faster on large inputs); `COAL_ENGINE=python|cpp` chooses.
- **Quieter, clearer engine output:** a wrong genetic code is reported once, before at most 20 per-codon warnings
  (`--max-warnings`, 0 = all; the `pal2nal` drop-in still prints them all); `--help` is grouped and short.
- **MCP results tell the truth:** `codon_align` returns `completed` (COAL finished), `all_ok` / `needs_attention`
  (whether every sequence matched its protein), `counts`, `suggested_codontable` and `next_step`, and a short summary
  of COAL's messages instead of raw stderr. The old `ok` field (which meant "finished") is gone. A missing file, an
  unknown genetic code or an unsupported format is now a tool error (`isError`), and the arguments have descriptions
  and allowed values. The server reports its version.
- **Windows, from source:** `desktop\run_from_source.bat` (the window) and `desktop\setup_mcp.bat` (the MCP server
  for Claude, printing ready-to-paste registrations) use separate environments, pick Python 3.12 when available,
  check the Python version and repair half-finished installs.
- **Windows app:** the MCP endpoint is off by default (it was on a random port no client could be registered with);
  `COAL.exe --mcp` turns it on.

## 0.4.1

The first release with a working Windows installer (the 0.4.0 build failed, so 0.4.0 has no downloads).

- **Clearer app:** a short header, section headings (Files, Genetic code, Output, Results), the genetic code first,
  a "Start alignment" button, and no empty results area before a run. When the genes fit another genetic code, the
  result says "Wrong genetic code?" with a one-click re-run instead of a page of warnings. Long lists of genes with the
  same problem are folded. One settings place (Gradio's own settings link is hidden). Amber accent colour.
- **Large inputs:** the alignment view is an honest preview at scale (problems first, overview strips, centre on any
  column), and a Settings tab remembers your defaults.
- **Better messages:** errors say what to do next (e.g. "Add --keep-going to write the other sequences"), swapped or
  unaligned input files are recognised, missing files are named plainly, and frameshift locations are only given as
  "at column N" when the evidence confirms it.
- **Windows app fixes:** the packaged app no longer crashes on start-up (UTF-8 output), and file paths with any
  characters work.
- **Examples** set their own genetic code, so one example can't leave the wrong code behind for the next.
- **Tested on real data:** 17 published codon alignments are rebuilt exactly; 12 ship as tests.

## 0.4.0

COAL turns a protein alignment plus the unaligned coding DNA into a codon alignment, like PAL2NAL v14, with the same
output apart from a short list of documented, deliberate differences (see the README).

- **Windows app:** an installer (or portable zip) with a local GUI: a plain-language health check for every gene, an
  alignment view that marks where a gene breaks, one-click re-run with a suggested genetic code, and a ready-to-paste
  methods sentence. Your data stays on your computer.
- **All four PAL2NAL output formats** (`-output fasta|clustal|paml|codon`), `-blockonly`, and a drop-in `pal2nal`
  command that accepts an unchanged PAL2NAL command line.
- **Run report** (`--report FILE.json` / `--json`) with each gene's status, the genetic-code suggestion and the methods
  sentence; `--keep-going` writes every other gene when one can't be aligned (exit 2); `--strict` refuses
  probably-wrong data; `--max-warnings`, `--progress` and `--resume` for very large runs.
- **Safer pairing:** IDs that match only after removing a version suffix, case or `sp|…|` fields are paired by ID,
  and COAL says so.
- **AI agents:** an MCP server with two tools, `codon_align` and `explain_report`.
- **C++ engine** (`coal-cpp`, built from source with CMake): same options and byte-identical output, 3-6x faster
  than the Python engine in our measurements.
- **Tested:** 285 stored checks against PAL2NAL v14's own output and 20,000 random comparisons with the real PAL2NAL
  v14, with 0 unexplained differences.

## 0.3.0

- **Large inputs:** a streaming pipeline with constant memory, automatic multi-process `-j` (only when each worker
  gets enough work) and an optional numpy fast path (`".[fast]"`); 500,000 sequences: 2.5 GB -> 143 MB of memory.
- All four PAL2NAL output formats and `-blockonly`, `--strict`, the drop-in `pal2nal` command, the GUI, MCP tools
  and the Windows app (see README.md).

## 0.2.0

- NCBI genetic codes 1-23 (`-codontable`), extra sequence around the gene (UTRs, stop codon) skipped, PAL2NAL's
  frameshift digits, `-nogap` and `-nomismatch`, and PAL2NAL's handling of codons that don't match.
- About twice as fast; on Windows 29x faster than PAL2NAL on 10,000 sequences and 119x on 10,000-codon genes,
  with identical output.

## 0.1.0

- First version: a pure-Python, dependency-free `coal` command that builds a codon alignment from a protein
  alignment and coding DNA, checked against stored PAL2NAL v14 outputs.
