# MCP demo: "which of my orthologs failed, and why?" (60 seconds)

A stage script for showing COAL as a tool an AI assistant can use. Everything below the script is the
real output of the tool functions on this machine (COAL engine 0.3.0), regenerated with
`python demo/make_mcp_demo.py`. Newer engines give the same story with more detail.

## Setup (before going on stage)

1. `pip install -e ".[mcp]"` in the COAL repo (or `.[gui]` to use the web GUI's MCP endpoint instead).
2. Add COAL to the assistant. Claude Code (run in the COAL folder, so `examples/...` paths resolve):

   ```
   claude mcp add --transport stdio coal -- python -m coal_app --mcp
   ```

   Claude Desktop (`claude_desktop_config.json`). Use the full path of the Python that has COAL installed
   (`python -c "import sys; print(sys.executable)"`), and give the assistant full paths to the example files:

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
| 0:20 | Tool call `explain_report` | "COAL found the problem gene, and it says *where*: a missing base at column 102." |
| 0:30 | Assistant's answer | "The assistant can explain it because COAL explains it: plain words, not a stack trace." |
| 0:38 | Second prompt (mito) | "Now mitochondrial genes, with the wrong settings, as people often do." |
| 0:45 | `codon_align` -> hint -> `codon_align` with code 2 | "COAL spots the wrong genetic code, suggests code 2, the assistant re-runs it. Everything fits." |
| 0:55 | Methods text in the answer | "And it hands you the methods sentence for your paper. That's COAL." |

**Prompt 1:**
> Make codon alignments for the orthologs in examples/one_broken_gene (protein.aln.fasta and dna.fasta)
> and tell me which failed and why.

**Prompt 2:**
> Same for examples/mito.

Expected assistant answer to prompt 1 (paraphrased): "OG0001_sp14 failed: one DNA base is missing at
alignment column 102 (a frameshift), so the codons after it are out of step. Check that gene's sequence
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
  "ok": false,
  "exit_code": 1,
  "output_path": null,
  "report_path": "<tmp>/protein.report.json",
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
  "methods_text": "",
  "coal_command": "coal protein.aln.fasta dna.fasta -o protein.codon.fasta",
  "pal2nal_command": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta > protein.codon.fasta"
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
  "ok": true,
  "exit_code": 0,
  "output_path": "<tmp>/protein.codon.fasta",
  "report_path": "<tmp>/protein.report.json",
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
      "text": "All sequences fit genetic code 2 (Vertebrate mitochondrial) without any mismatches, compared with 0 under code 1. Re-run with code 2?",
      "suggested_codontable": 2
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
    "... 21 more"
  ],
  "methods_text": "Codon alignments were built with COAL v0.3.0 (https://github.com/chloe-newt/COAL), which places each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method of PAL2NAL; Suyama et al. 2006), using the standard genetic code (NCBI translation table 1). Of 24 sequences, 0 matched their protein exactly, 24 contained codons that did not match the aligned amino acid (kept in the alignment). Command: coal protein.aln.fasta dna.fasta -o protein.codon.fasta",
  "coal_command": "coal protein.aln.fasta dna.fasta -o protein.codon.fasta",
  "pal2nal_command": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta > protein.codon.fasta"
}
```

### 4. `codon_align` on mito again with `codontable=2` (from the hint)

```json
{
  "ok": true,
  "exit_code": 0,
  "output_path": "<tmp>/protein.codon.fasta",
  "report_path": "<tmp>/protein.report.json",
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
  "methods_text": "Codon alignments were built with COAL v0.3.0 (https://github.com/chloe-newt/COAL), which places each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method of PAL2NAL; Suyama et al. 2006), using the vertebrate mitochondrial genetic code (NCBI translation table 2). Of 24 sequences, 24 matched their protein exactly. Command: coal protein.aln.fasta dna.fasta -codontable 2 -o protein.codon.fasta",
  "coal_command": "coal protein.aln.fasta dna.fasta -codontable 2 -o protein.codon.fasta",
  "pal2nal_command": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta -codontable 2 > protein.codon.fasta"
}
```
