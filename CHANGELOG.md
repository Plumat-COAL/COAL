# Changelog

## Unreleased

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
  are written through, and FIFOs, devices and hard-linked files are written in place.
- **A killed worker no longer hangs the Python engine** (robustness review): when the out-of-memory killer (or
  `kill -9`) took one worker of a `-j` run, the main process could block forever in the pool's shutdown. It now
  ends the pool without waiting and stops with `coal: error: a worker process died (out of memory?)` (exit 1, report
  code `worker_died`, no partial output, no temp files). Workers also let SIGTERM end them (the inherited handler
  used to be swallowed inside the work) and **die with the main process** (`PR_SET_PDEATHSIG` on Linux, a parent
  watchdog everywhere), so a `kill -9` of `coal` no longer leaves workers holding its stdout open and a pipeline
  reading it waiting for EOF forever.
- **Checkpoints are checked before `--resume` trusts them** (robustness review): an edited or corrupted
  `<out>.coal-checkpoint.json` (`completed` beyond the input or negative, a string or float where an integer belongs,
  a list at the top level, a missing key, a deleted entries file, a partial file shorter than recorded, a truncated
  or deeply nested file, a folder at the checkpoint's path) made both engines write a truncated output with exit 0,
  or crash with a traceback / raw C++ text. Both engines now check every field's type and range and the partial
  files' sizes, treat any problem as the existing "unreadable checkpoint" notice and start a normal run; a
  checkpoint that cannot be written (`cannot write the checkpoint ...: is a directory`) is a notice and the run
  goes on without one; and `coal-cpp`'s JSON parser has a nesting limit (64) instead of overflowing the stack.
- **A checkpoint written by the other engine** (robustness review): `coal --resume` now picks the engine that wrote
  the checkpoint (so a `coal-cpp` installed or removed between the run and the resume does not get in the way),
  and when an engine is asked to resume the other one's run it says so instead of blaming the inputs: `out.fa was
  interrupted by the Python engine (coal-python), which must resume it: run again with COAL_ENGINE=python, or start
  a new run without --resume`.
- **`-o /dev/fd/1` and `-o /dev/stdout` work when stdout is a pipe** (robustness review): both engines resolved
  the path to `pipe:[N]` and failed with `cannot write /dev/fd/1: No such file or directory`. An existing target
  that is not a regular file is now written through before any path resolution.
- **`--json | head` is a silent exit 141** (robustness review) in both engines (the Python engine printed a
  traceback and exited 1, `coal-cpp` exited 0), also for `--batch ... --json`.
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
  `os.PathLike` now means that one file (the string used to be iterated character by character).
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
- **Documented for 1.0:** `docs/python-api.md` (the public Python API: `coal.codon_align`, the new re-export
  `coal.align_files` = the `coal` command as a function, `coal.CoalError` with its `code`/`seq_id`/`index`, and the
  semantic-versioning promise for command-line flags, output bytes, exit codes, the report schema and these
  functions), `docs/report-schema.md` with a JSON Schema (`docs/report.schema.json`, draft 2020-12) for the run report
  `coal.report/1` (every field, status, pairing mode, message and hint code; additive changes only within `/1`),
  `CONTRIBUTING.md`, GitHub issue forms and a pull-request template. New tests run the documented examples
  (`tests/test_python_api.py`) and validate freshly generated reports from both engines against the schema
  (`tests/test_report_schema.py`); `tests/report_example.json` is now a real `examples/one_broken_gene --keep-going`
  report. No change to COAL's behaviour or output.
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
- **Batch mode (`coal --batch JOBS`, #15):** thousands of small conversions in one process, one `coal` command line
  per line of JOBS (`-` = stdin), the command line's options as the defaults of every job, `--batch-out-dir DIR` for
  jobs without `-o`, a `coal.batch/1` report with `--report`/`--json`, per-job stderr prefixes
  (`coal: [k/n PEP] warning: ...`) and a `coal: batch:` summary line. A failing job does not stop the batch (exit 1
  at the end; 2 if jobs only left sequences out; mistakes in JOBS are reported with their line number before
  anything runs). Every output is byte-identical to the single run; about 40× less wall time on tiny jobs. Python
  only: `coal-cpp --batch` is a clear error and batches are never dispatched to it. `--resume` and cross-job
  parallelism are not supported in a batch yet.
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
