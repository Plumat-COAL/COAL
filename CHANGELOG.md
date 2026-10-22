# Changelog

## Unreleased

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
