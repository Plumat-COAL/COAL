---
name: coal
description: Build codon alignments with COAL (a fast, MIT-licensed PAL2NAL replacement) from a protein multiple sequence alignment plus the unaligned coding DNA of the same genes. Use when someone needs a codon alignment, back-translated alignment, "pal2nal", input for dN/dS (PAML codeml, HyPhy), or asks why protein and DNA "don't match".
---

<!-- CHANGED(APP-C10, APP-C11): was documenting upcoming engine options as if present, "exit 2 = partial
success" and both extras in one environment; now current options only, upcoming ones in their own section,
exit 2 = usage error today, extras in separate environments. See overnight/app/CHANGES.md -->

# COAL: codon alignment from a protein alignment + coding DNA

## When to use it

- You have **aligned proteins** (from MAFFT, MUSCLE, Clustal...) and the **unaligned CDS/DNA** of the same genes,
  and need the DNA aligned codon by codon (the input for dN/dS, selection scans, codon-aware trees).
- Someone mentions PAL2NAL: COAL accepts the same options (`-codontable`, `-nogap`, `-nomismatch`,
  `-output fasta`) and gave the same FASTA output as PAL2NAL v14 in all 147 byte-for-byte comparisons of its
  test suite; it is much faster on big inputs. <!-- CHANGED(R2-D2): was 143 -->
- Not for: aligning sequences from scratch (run an aligner on the proteins first), non-coding DNA, or RNA genes.

## Install

```
pip install -e .            # from the COAL repo: pure Python 3.10+, no dependencies
pip install -e .[gui]       # optional: web GUI, which also serves the MCP tools (Gradio)
pip install -e .[mcp]       # optional: stand-alone stdio MCP server only
```

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
```

<!-- CHANGED(R2-A10, R2-D5, R2-D7, R2-E4): was "-q silences them", "unrelated sequences can still succeed",
"always count the warnings per gene" (none are printed with -nomismatch) and exit codes 0/1/2 only.
See overnight/integration/CHANGES.md -->
Output is FASTA only (`-output fasta`, the default). Per-codon warnings (`coal: warning: <id> pepAlnPos N: ...`)
go to stderr, and so do **notices about probably-wrong data** (below). `-q` silences warnings and notices alike
(errors are always printed), so **don't use `-q` when you need to know about bad genes**; use `--strict` to stop on
them instead. <!-- CHANGED(DEC-2): was "-q silences those, but notices ... are always printed" --> The notices:
`<id>: no 10-residue segment of the protein was found in the DNA; codons were assigned by position (K of N
checkable codons do not match; A contain N or other ambiguous bases)` (the last part only when A > 0) and `only K of N protein IDs are DNA IDs; pairing sequences by order instead (e.g. protein 'A'
gets DNA 'B')` (some proteins are then almost certainly paired with the wrong DNA: fix the IDs).
<!-- CHANGED(DEC-1): was "that is an error" for the next two cases; now PAL2NAL's output plus a notice, an error only
with --strict. See overnight/integration/CHANGES.md -->
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
COAL's messages).

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
    (paired with another gene's DNA by order). Don't use such a row. <!-- CHANGED(DEC-1) -->
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
| `number of sequences differs` | Protein and DNA files list different genes | Make both files contain the same genes |
| `protein and DNA are inconsistent` for the **first** gene, or for many genes | IDs don't match (e.g. `AT1G01010.1` vs `AT1G01010`), so COAL paired by order | Make IDs identical (strip `.1` version suffixes) |
| `protein and DNA are inconsistent` for **one** gene | Frameshift (missing/extra base), wrong transcript, or truncated DNA | Check that gene's DNA; drop it and re-run (the current engine stops at the first failing gene) |
| `invalid protein alignment characters` | Non-amino-acid symbols in the alignment | Remove them; `.` and `-` are both gaps |
| Stop codons in output columns | Internal stops or terminal stop kept | `-nogap` removes stop columns |

Never "fix" a mismatch by editing sequences silently: report which genes failed and why.

## After the engine update (only if `python -m coal --help` lists these)

These are being built on the engine branch; do not use them unless `--help` shows them.

```
coal pep.aln nuc.fasta -output paml -o codon.paml        # also -output clustal / codon
coal pep.aln nuc.fasta --keep-going --report run.json -o codon.fasta   # pipelines: keep going, write a report
```

- With `--keep-going`, exit code 2 also means **partial success**: some genes failed, the rest were written
  (check the report). Without `--keep-going` support, exit 2 is only ever a usage error.
- `--report run.json` writes the `coal.report/1` report directly.
- ID normalisation: `pairing: normalised` (IDs matched after tidying, e.g. `.1` removed).
