---
name: coal
description: Build codon alignments with COAL (a fast, MIT-licensed PAL2NAL replacement) from a protein multiple sequence alignment plus the unaligned coding DNA of the same genes. Use when someone needs a codon alignment, back-translated alignment, "pal2nal", input for dN/dS (PAML codeml, HyPhy), or asks why protein and DNA "don't match".
---

# COAL: codon alignment from a protein alignment + coding DNA

## When to use it

- You have **aligned proteins** (from MAFFT, MUSCLE, Clustal...) and the **unaligned CDS/DNA** of the same genes,
  and need the DNA aligned codon by codon (the input for dN/dS, selection scans, codon-aware trees).
- Someone mentions PAL2NAL: COAL accepts the same options (`-codontable`, `-nogap`, `-nomismatch`, `-output`)
  and gives byte-identical FASTA output, much faster on big inputs.
- Not for: aligning sequences from scratch (run an aligner on the proteins first), non-coding DNA, or RNA genes.

## Install

```
pip install -e .            # from the COAL repo: pure Python 3.10+, no dependencies
pip install -e .[gui]       # optional: web GUI + MCP server (Gradio)
pip install -e .[mcp]       # optional: stdio MCP server only
```

Check: `python -m coal --version`. Use `python -m coal` if the `coal` command is not on PATH (common on Windows).

## The 5 most common commands

```
coal pep.aln nuc.fasta -o codon.fasta                    # 1. basic run (nuclear genes, standard code)
coal pep.aln nuc.fasta -codontable 2 -o codon.fasta      # 2. vertebrate mitochondrial genes
coal pep.aln nuc.fasta -nogap -o codon.fasta             # 3. drop gap/stop columns (for codeml)
coal pep.aln nuc.fasta -output paml -o codon.paml        # 4. PAML format (if `coal --help` lists it)
coal pep.aln nuc.fasta --keep-going --report run.json -o codon.fasta   # 5. pipelines: keep going, write a report
```

Options 4 and 5 are in newer engines: check `python -m coal --help` first. Without them, use FASTA output and
read the warnings on stderr. Exit codes: 0 = success, 1 = error (nothing written), 2 = partial success
(with `--keep-going`: some genes failed, the rest were written).

Genetic codes (`-codontable N`): 1 standard (nuclear), 2 vertebrate mitochondrial, 3 yeast mito,
4 mold/protozoan mito and Mycoplasma, 5 invertebrate mito, 9 echinoderm/flatworm mito, 11 bacteria, archaea
and plastids; also 6, 10, 12-16, 21-23.

GUI for people: `python -m coal_app` (opens on http://127.0.0.1:7860). AI agents: `python -m coal_app --mcp`
gives two MCP tools, `codon_align(...)` and `explain_report(report_path)`.

## How to read the report (`--report run.json`, schema `coal.report/1`)

- `summary`: `total`, `ok`, `warning`, `failed`, and `pairing` = how proteins met their DNA:
  `id` (exact IDs), `normalised` (IDs matched after tidying, e.g. `.1` removed), `order` (by position;
  risky, check it).
- `sequences[]`: one entry per sequence with `status` `ok` / `warning` / `failed` and `messages[]`
  (`code`, `text`, `column` = protein-alignment column, 1-based).
  - `warning` + `mismatch`: a codon doesn't code for its amino acid; the row is still in the output.
  - `failed`: COAL could not fit the protein onto the DNA at all; that row is not in the output.
- `hints[]`: suggestions. `genetic_code` comes with `suggested_codontable`; re-run with it.
- `output`: path, format and sha256 of the alignment; `null` if nothing was written.
- `command` / `pal2nal_equivalent`: the exact command lines, for methods sections.

For a plain-English summary call the MCP tool `explain_report(report_path)`.

## Common failures and fixes

| What you see | Likely cause | Fix |
|---|---|---|
| Many `X does not correspond to TGA` warnings (W vs TGA), or M vs ATA | Mitochondrial genes read with the standard code | `-codontable 2` (vertebrates), 5 (invertebrates), 4, 9...; follow the `genetic_code` hint |
| `number of sequences differs` | Protein and DNA files list different genes | Make both files contain the same genes |
| `protein and DNA are inconsistent` for the **first** gene, or for many genes | IDs don't match (e.g. `AT1G01010.1` vs `AT1G01010`), so COAL paired by order | Make IDs identical (strip `.1` version suffixes), or use an engine with ID normalisation (`pairing: normalised`) |
| `protein and DNA are inconsistent` for **one** gene | Frameshift (missing/extra base), wrong transcript, or truncated DNA | Check that gene's DNA; drop it, or use `--keep-going` to align the rest |
| `invalid protein alignment characters` | Non-amino-acid symbols in the alignment | Remove them; `.` and `-` are both gaps |
| Stop codons in output columns | Internal stops or terminal stop kept | `-nogap` removes stop columns |

Never "fix" a mismatch by editing sequences silently: report which genes failed and why.
