# MCP demo: "which of my orthologs failed, and why?" (60 seconds)

A stage script for showing COAL as a tool an AI assistant can use. Everything below the script is the
real output of the tool functions, regenerated with `python demo/make_mcp_demo.py`.
<!-- CHANGED(IMP-M-4): version-free (was "COAL engine 0.3.0"); setup points to desktop/setup_mcp.bat -->

## Setup (before going on stage)

1. `pip install -e ".[mcp]"` in the COAL repo, in its own virtual environment (not the one with `".[gui]"`).
   On Windows, `desktop\setup_mcp.bat` does this and prints the two registrations below, filled in.
2. Add COAL to the assistant. Claude Code (use the full path of that environment's Python if it is not the
   `python` on your PATH):

   ```
   claude mcp add --transport stdio coal -- python -m coal_app --mcp
   ```

   Claude Desktop (`claude_desktop_config.json`). Use the full path of the Python that has COAL installed
   (`python -c "import sys; print(sys.executable)"`). Give the assistant full paths to the files in both cases:

   ```json
   {
     "mcpServers": {
       "coal": {
         "command": "/full/path/to/python",
         "args": ["-m", "coal_app", "--mcp"]
       }
     }
   }
   ```

   With the web GUI running (`python -m coal_app`), MCP clients that speak SSE can use
   `http://127.0.0.1:7860/gradio_api/mcp/sse` instead.
3. Check the two tools are listed (`codon_align`, `explain_report`). Have `examples/` open in a file browser.

## The script

| Time | On screen | Say |
|---|---|---|
| 0:00 | Chat window | "Bioinformatics tools are built for people at a terminal. COAL is also built for AI assistants." |
| 0:05 | Type the prompt below | "I'll ask in plain English, like a biologist would." |
| 0:12 | Tool call `codon_align` on one_broken_gene | "It calls COAL directly. No copy-pasting commands, no file-format guessing." |
| 0:20 | Tool call `explain_report` | "COAL found the problem gene, and COAL's app estimates *where*: probably a missing base at column 102." |
| 0:30 | Assistant's answer | "The assistant can explain it because COAL explains it: plain words, not a stack trace." |
| 0:38 | Second prompt (mito) | "Now mitochondrial genes, with the wrong settings, as people often do." |
| 0:45 | `codon_align` -> hint -> `codon_align` with code 2 | "COAL spots the wrong genetic code, suggests code 2, the assistant re-runs it. Everything fits." |
| 0:55 | Methods text in the answer | "And it hands you the methods sentence for your paper. That's COAL." |

**Prompt 1:**
> Make codon alignments for the orthologs in examples/one_broken_gene (protein.aln.fasta and dna.fasta)
> and tell me which failed and why.

**Prompt 2:**
> Same for examples/mito.

Expected assistant answer to prompt 1 (paraphrased): "OG0001_sp14 failed: COAL's app estimates that one DNA
base is missing at alignment column 102 (a frameshift), so the codons after it are out of step. Check that gene's sequence
(likely a sequencing or assembly error). COAL stopped there, so no alignment was written yet: drop that gene
and re-run, or (with a COAL that has `--keep-going`) align the other 19 now." With a `--keep-going` engine
the first call already returns the alignment of the 19 good genes, and the answer starts "19 of 20 are fine".

## Real tool outputs

### 1. `codon_align` on one_broken_gene

Call:
```json
{
  "protein_alignment_path": "<COAL>/examples/one_broken_gene/protein.aln.fasta",
  "dna_fasta_path": "<COAL>/examples/one_broken_gene/dna.fasta"
}
```
Result:
```json
{
  "completed": false,
  "all_ok": false,
  "needs_attention": true,
  "counts": {
    "total": 20,
    "ok": 0,
    "warning": 0,
    "failed": 1,
    "not_checked": 19
  },
  "suggested_codontable": null,
  "next_step": "COAL stopped before writing an alignment: call explain_report(report_path) for the cause and fix.",
  "output_path": null,
  "report_path": "<tmp>/protein.report.json",
  "exit_code": 1,
  "summary": {
    "total": 20,
    "ok": 0,
    "warning": 0,
    "failed": 1,
    "pairing": "id",
    "not_checked": 19
  },
  "summary_text": "20 sequences: 0 OK, 0 with warnings, 1 failed, 19 not checked.",
  "hints": [],
  "problems": [
    {
      "id": "OG0001_sp14",
      "status": "failed",
      "explanation": "Frameshift at column 102: one base is missing in the DNA, so from here on the codons are read out of step. Check the DNA near there (a sequencing or assembly error?)."
    }
  ],
  "problems_total": 1,
  "engine_messages": [
    "coal: error: OG0001_sp14: protein and DNA are inconsistent (no match found)"
  ],
  "warning_lines": 0,
  "methods_text": "",
  "coal_command": "coal protein.aln.fasta dna.fasta -o protein.codon.fasta",
  "pal2nal_command": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta > protein.codon.fasta",
  "pal2nal_caveats": []
}
```

### 2. `explain_report` on that report

```
20 sequences: 0 OK, 0 with warnings, 1 failed, 19 not checked.
Proteins and DNA were paired by matching their IDs exactly.
No alignment was written, because COAL stopped at an error.

Failed sequences:
- OG0001_sp14: Frameshift at column 102: one base is missing in the DNA, so from here on the codons are read out of step. Check the DNA near there (a sequencing or assembly error?).

How to fix a failed sequence: check it is the right DNA for that protein (same gene, same species), look for a missing or extra base near where it stops matching, and make sure the genetic code is right for where the gene comes from (for example code 2 for vertebrate mitochondria).
```

### 3. `codon_align` on mito (default genetic code 1)

```json
{
  "completed": true,
  "all_ok": false,
  "needs_attention": true,
  "counts": {
    "total": 24,
    "ok": 0,
    "warning": 24,
    "failed": 0,
    "not_checked": 0
  },
  "suggested_codontable": 2,
  "next_step": "Re-run codon_align with codontable=2 (see hints), then check again.",
  "output_path": "<tmp>/protein.codon.fasta",
  "report_path": "<tmp>/protein.report.json",
  "exit_code": 0,
  "summary": {
    "total": 24,
    "ok": 0,
    "warning": 24,
    "failed": 0,
    "pairing": "id"
  },
  "summary_text": "24 sequences: 0 OK, 24 with warnings, 0 failed.",
  "hints": [
    {
      "code": "genetic_code",
      "text": "Under genetic code 2 (Vertebrate mitochondrial), all 24 sequences match their protein; under code 1, none do. Codes 5 and 13 fit these sequences too: choose the code of the organism and organelle the genes come from. Re-run with code 2?",
      "suggested_codontable": 2,
      "checked": 24,
      "fit": 24,
      "fit_current": 0
    }
  ],
  "problems": [
    {
      "id": "Homo_COX1",
      "status": "warning",
      "explanation": "7 codons do not match their amino acid (first at column 1): the protein says M but the DNA has ATA. Kept in the alignment; check the genetic code."
    },
    {
      "id": "Pan_COX1",
      "status": "warning",
      "explanation": "4 codons do not match their amino acid (first at column 1): the protein says M but the DNA has ATA. Kept in the alignment; check the genetic code."
    },
    {
      "id": "Mus_COX1",
      "status": "warning",
      "explanation": "13 codons do not match their amino acid (first at column 1): the protein says M but the DNA has ATA. Kept in the alignment; check the genetic code."
    },
    "... 7 more"
  ],
  "problems_total": 24,
  "engine_messages": [],
  "warning_lines": 179,
  "methods_text": "Codon alignments were built with COAL v0.3.0 (https://github.com/plumat-coal/COAL), which places each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method of PAL2NAL; Suyama et al. 2006), using the standard genetic code (NCBI translation table 1). Of 24 sequences, 0 matched their protein exactly, and 24 contained codons that did not match the aligned amino acid (kept in the alignment). Command: coal protein.aln.fasta dna.fasta -o protein.codon.fasta",
  "coal_command": "coal protein.aln.fasta dna.fasta -o protein.codon.fasta",
  "pal2nal_command": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta > protein.codon.fasta",
  "pal2nal_caveats": []
}
```

### 4. `codon_align` on mito again with `codontable=2` (from the hint)

```json
{
  "completed": true,
  "all_ok": true,
  "needs_attention": false,
  "counts": {
    "total": 24,
    "ok": 24,
    "warning": 0,
    "failed": 0,
    "not_checked": 0
  },
  "suggested_codontable": null,
  "next_step": "All sequences match their protein; the alignment at output_path is ready to use.",
  "output_path": "<tmp>/protein.codon.fasta",
  "report_path": "<tmp>/protein.report.json",
  "exit_code": 0,
  "summary": {
    "total": 24,
    "ok": 24,
    "warning": 0,
    "failed": 0,
    "pairing": "id"
  },
  "summary_text": "24 sequences: 24 OK, 0 with warnings, 0 failed.",
  "hints": [],
  "problems": [],
  "problems_total": 0,
  "engine_messages": [],
  "warning_lines": 0,
  "methods_text": "Codon alignments were built with COAL v0.3.0 (https://github.com/plumat-coal/COAL), which places each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method of PAL2NAL; Suyama et al. 2006), using the vertebrate mitochondrial genetic code (NCBI translation table 2). Of 24 sequences, 24 matched their protein exactly. Command: coal protein.aln.fasta dna.fasta -codontable 2 -o protein.codon.fasta",
  "coal_command": "coal protein.aln.fasta dna.fasta -codontable 2 -o protein.codon.fasta",
  "pal2nal_command": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta -codontable 2 > protein.codon.fasta",
  "pal2nal_caveats": []
}
```
