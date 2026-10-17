# Changelog

## Unreleased

<!-- CHANGED(IMP-M-4): Unreleased, 0.2.0 and 0.1.0 added from the git history -->
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
