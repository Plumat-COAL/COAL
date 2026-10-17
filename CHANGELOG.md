# Changelog

## 0.4.0

<!-- CHANGED(MRG-8): new file (scripts/release_notes.py takes the release notes from this section). -->

- Merged the extras and C++ branch (overnight/core) into the reviewed engine:
  - `--report FILE.json` / `--json`: the run report (schema `coal.report/1`) written by the engine itself, with the
    same per-gene statuses, messages, genetic-code suggestion and methods sentence as the GUI's report;
    `coal: summary:` lines at the end of a run with problems.
  - `--keep-going` (skip genes that cannot be converted, exit 2), `--max-warnings N`, `--progress`, `--resume`
    (checkpoints; resumed output byte-identical).
  - Pairing by normalised IDs (version suffix, case, `sp|…|` fields) where PAL2NAL's file order would give a protein
    the DNA of another gene (a documented deliberate difference); unchanged where the order is right.
  - A new reference case, `frameshift_codon_line` (285 stored PAL2NAL checks).
  - The C++ engine (`cpp/`, CMake) and its tests; it is being brought up to the Python engine's behaviour.
- Version 0.4.0.

## 0.3.0

- All four PAL2NAL output formats and `-blockonly`, `--strict`, the drop-in `pal2nal` command, the GUI, MCP tools
  and the Windows app (see README.md).
