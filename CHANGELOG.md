# Changelog

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

<!-- CHANGED(MRG-8): new file (scripts/release_notes.py takes the release notes from this section). -->
<!-- CHANGED(SHIP-6): was the internal merge log ("Merged the extras and C++ branch ..."); now written for users -->

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

- All four PAL2NAL output formats and `-blockonly`, `--strict`, the drop-in `pal2nal` command, the GUI, MCP tools
  and the Windows app (see README.md).
