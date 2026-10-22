---
name: coal
description: Build codon alignments with COAL (a fast, MIT-licensed PAL2NAL replacement) from a protein multiple sequence alignment plus the unaligned coding DNA of the same genes. Use when someone needs a codon alignment, back-translated alignment, "pal2nal", input for dN/dS (PAML codeml, HyPhy), or asks why protein and DNA "don't match".
---

# COAL: codon alignment from a protein alignment + coding DNA

## When to use it

- You have **aligned proteins** (from MAFFT, MUSCLE, Clustal...) and the **unaligned CDS/DNA** of the same genes,
  and need the DNA aligned codon by codon (the input for dN/dS, selection scans, codon-aware trees).
- Someone mentions PAL2NAL: COAL accepts the same options (`-codontable`, `-nogap`, `-nomismatch`, `-blockonly`,
  `-output clustal|paml|fasta|codon`) and gave the same output as PAL2NAL v14 in all 238 byte-for-byte comparisons
  of its test suite (all four formats), plus a random differential test against the real PAL2NAL; it is much faster
  on big inputs. The `pal2nal` command takes an unchanged PAL2NAL command line (default clustal, like PAL2NAL).
- Not for: aligning sequences from scratch (run an aligner on the proteins first), non-coding DNA, or RNA genes.

## Install

```
pip install -e .                              # from the COAL repo: the engine, pure Python 3.10+, no dependencies
pip install -e . -e "apps/coal-app[gui]"      # optional: web GUI, which also serves the MCP tools (Gradio)
pip install -e . -e "apps/coal-app[mcp]"      # optional: stand-alone stdio MCP server only
```

The engine (`coal-align`) is all a pipeline needs; the GUI and the MCP server are the second package,
`coal-align-app` (`apps/coal-app`). Without a checkout:
`pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz"`
(add `"coal-align-app[mcp] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"` for
the MCP server).

Install `gui` and `mcp` in **separate environments** (virtualenvs): Gradio's own MCP support may need an
older `mcp` than the `mcp` extra pins. Check: `python -m coal --version`. Use `python -m coal` if the `coal`
command is not on PATH (common on Windows).

## First: check what this engine supports

Run `python -m coal --help` before using any option. Engines differ (see "After the engine update" below);
an option the engine doesn't know is a **usage error**: argparse prints it and exits with code **2**.

## The most common commands (current engine)

```
coal pep.aln nuc.fasta -o codon.fasta                    # 1. basic run (nuclear genes, standard code)
coal pep.aln nuc.fasta -codontable 2 -o codon.fasta      # 2. vertebrate mitochondrial genes
coal pep.aln nuc.fasta -nogap -o codon.fasta             # 3. drop gap/stop columns (for codeml)
coal pep.aln nuc.fasta -nomismatch -o codon.fasta        # 4. drop columns where a codon doesn't fit
coal pep.aln nuc.fasta -j 4 -o codon.fasta               # 5. big inputs: up to 4 worker processes
coal pep.aln nuc.fasta -output paml -o codon.paml        # 6. PAML/codeml input (also -output clustal / codon)
coal pep.aln nuc.fasta --strict -o codon.fasta           # 7. pipelines: stop on probably-wrong data
```

Output formats: `-output fasta` (the default of `coal`), `clustal`, `paml`, `codon`, as PAL2NAL v14 writes them.
COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate differences (newer NCBI genetic-code tables for codes 3, 10, 11 and 13; `-nogap` removing the selected code's stop codons; the selected code's start codons for the first Met; pairing by IDs that match only after normalising where PAL2NAL v14's file order would give a protein another gene's DNA; and inputs where PAL2NAL v14 writes malformed output or fails). Tested: all 285 stored reference checks pass (238 byte-identical, 4 of them with deliberately different warnings; 24 where both tools fail; 23 deliberate output or exit-status differences: 6 deliberate errors, 3 different outputs, 14 `-nogap` checks with a non-standard code), and in 20,000 random comparisons with the real PAL2NAL v14 (5,000 cases x 4 output formats) 13,826 were byte-identical, 4,681 both failed and 1,493 were deliberate differences, each checked to be exactly the documented change; 0 unexplained differences. (List in docs/differences.md.) (`codon` cannot be combined with `-nogap`/`-nomismatch`/`-blockonly`, as in PAL2NAL). `-blockonly` keeps
only the columns marked `#` under the blocks of a CLUSTAL input. Per-codon warnings (`coal: warning: <id> pepAlnPos N: ...`)
go to stderr, and so do **notices about probably-wrong data** (below). `-q` silences warnings and notices alike
(errors are always printed), so **don't use `-q` when you need to know about bad genes**; use `--strict` to stop on
them instead. The notices:
`<id>: no 10-residue segment of the protein was found in the DNA; codons were assigned by position (K of N
checkable codons do not match; A contain N or other ambiguous bases)` (the last part only when A > 0) and `only K of N protein IDs are DNA IDs; pairing sequences by order instead (e.g. protein 'A'
gets DNA 'B')` (some proteins are then almost certainly paired with the wrong DNA: fix the IDs).
Like PAL2NAL, the engine writes an alignment row even for probably-wrong data, and then prints a notice:
`<id>: protein and DNA are probably unrelated: K of N checkable codons do not match the protein. ...` (more than
half of the codons made only of A/C/G/T don't match: wrong DNA or wrong genetic code) and `<id>: position of the
gene in the DNA is a guess: ... only K of N codons have no ambiguous bases, too few to tell where the gene starts
in the longer DNA. ...` (heavily masked gene plus extra DNA: the row may be shifted; trim the UTR). **Treat such a
gene's row as wrong** (the report marks it `failed`). Add **`--strict`** to make these cases, and a partial ID
match, errors instead (exit 1, nothing written); check `python -m coal --help` for it first.

With `-nomismatch` the engine, like PAL2NAL, prints **no** mismatch warnings, so stderr cannot tell you
which genes didn't fit. Run once without `-nomismatch` to see them, or use `codon_align`: it re-checks
without `-nomismatch` and reports per-gene warnings and `summary.columns_removed`.

Exit codes: 0 = success (warnings possible), 1 = error (nothing written; the first error stops the run),
2 = usage error (bad option or a file name starting with `-`: put `--` before the file names), 130 =
interrupted (Ctrl-C), 141 = output pipe closed (e.g. `| head`), 143 = stopped by SIGTERM, 129 = SIGHUP;
in all of these nothing is written. On Windows (expected, not yet verified there): Ctrl-C gives 130 and a
closed pipe 141, but a process killed from outside (Task Manager, `taskkill /F`) exits 1 without cleaning up.

Genetic codes (`-codontable N`): 1 standard (nuclear), 2 vertebrate mitochondrial, 3 yeast mito,
4 mold/protozoan mito and Mycoplasma, 5 invertebrate mito, 9 echinoderm/flatworm mito, 11 bacteria, archaea
and plastids; also 6, 10, 12-16, 21-23.

GUI for people: `python -m coal_app` (opens on http://127.0.0.1:7860 and serves MCP there too). AI agents:
`python -m coal_app --mcp` gives two MCP tools, `codon_align(...)` and `explain_report(report_path)`;
`codon_align` returns a `report_path` even with the current engine (the app builds the report from
COAL's messages). Give both tools **absolute paths**. In the result, `completed` only means COAL finished;
use the alignment only when `all_ok` is true. If `needs_attention` is true, read `counts`
(ok / warning / failed), `next_step` and `problems`, or call `explain_report(report_path)`; if
`suggested_codontable` is set, re-run `codon_align` with it. A missing file, an unknown genetic code or an
unsupported format is a tool error (isError), not a result.

## How to read the report (`report_path` from `codon_align`, schema `coal.report/1`)

- `summary`: `total`, `ok`, `warning`, `failed`, `not_checked` (genes COAL never reached because it stopped
  at an earlier error), and `pairing` = how proteins met their DNA: `id` (exact IDs) or `order` (by position,
  because not every protein ID is a DNA ID; risky, check it); `null` if COAL stopped before pairing.
- `sequences[]`: one entry per sequence with `status` `ok` / `warning` / `failed` and `messages[]`
  (`code`, `text`, `column` = protein-alignment column, 1-based).
  - `warning` + `mismatch`: a codon doesn't code for its amino acid; the row is still in the output.
  - `failed`: COAL could not fit the protein onto the DNA at all (that row is not in the output), or the row was
    written as PAL2NAL writes it but the data is probably wrong: first message `suspect_unrelated` (most codons
    don't match), `suspect_position` (start of the gene can't be told), or `suspect_pairing` with mismatches
    (paired with another gene's DNA by order). Don't use such a row.
- `summary.columns_removed` (only with `-nomismatch`): how many codon columns were removed because a codon
  didn't match; the genes with such codons are `warning`, not `ok`.
- `hints[]`: suggestions. `genetic_code` comes with `suggested_codontable`; re-run with it.
  `partial_id_match`: only some protein IDs were found in the DNA, so pairing was by order and is probably
  wrong for some genes: tell the user, don't trust the alignment. `engine_warning`: any other COAL warning.
- `output`: path, format and sha256 of the alignment; `null` if nothing was written.
- `command` / `pal2nal_equivalent`: the exact command lines, for methods sections.

For a plain-English summary call the MCP tool `explain_report(report_path)`.

## Common failures and fixes

| What you see | Likely cause | Fix |
|---|---|---|
| Many `X does not correspond to TGA` warnings (W vs TGA), or M vs ATA | Mitochondrial genes read with the standard code | `-codontable 2` (vertebrates), 5 (invertebrates), 4, 9...; follow the `genetic_code` hint |
| `number of sequences differs` | Protein and DNA files list different genes | Make both files contain the same genes (the message names an extra or missing record when the IDs show it) |
| `... looks like DNA, not a protein alignment` | The two files were given the other way round | Protein alignment first, DNA second: `coal PROTEIN_ALIGNMENT DNA_FASTA` |
| `protein alignment rows have different lengths ...: is X an aligned protein file` | The protein file is not aligned | Align it first (MAFFT, MUSCLE, ...) |
| `cannot open 'x.fa': no such file` | Wrong path | Check the file name and folder |
| `protein and DNA are inconsistent` for the **first** gene, or for many genes | IDs don't match (e.g. `AT1G01010.1` vs `AT1G01010`), so COAL paired by order | Make IDs identical (strip `.1` version suffixes) |
| `protein and DNA are inconsistent` for **one** gene | Frameshift (missing/extra base), wrong transcript, or truncated DNA | Check that gene's DNA (the summary line says where it breaks); drop it and re-run (the current engine stops at the first failing gene) |
| `invalid protein alignment characters` | Non-amino-acid symbols in the alignment | Remove them; `.` and `-` are both gaps |
| Stop codons in output columns | Internal stops or terminal stop kept | `-nogap` removes stop columns |

Never "fix" a mismatch by editing sequences silently: report which genes failed and why.

## After the engine update (only if `python -m coal --help` lists these)

These are being built on the engine branch; do not use them unless `--help` shows them.

```
coal pep.aln nuc.fasta --keep-going --report run.json -o codon.fasta   # pipelines: keep going, write a report
```

- With `--keep-going`, exit code 2 also means **partial success**: some genes failed, the rest were written
  (check the report). Without `--keep-going` support, exit 2 is only ever a usage error.
- `--report run.json` writes the `coal.report/1` report directly.
- ID normalisation: `pairing: normalised` (IDs matched after tidying, e.g. `.1` removed).
