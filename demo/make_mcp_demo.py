"""Regenerate demo/mcp_demo.md by really calling the MCP tool functions on the examples.

    python demo/make_mcp_demo.py

Temporary paths in the outputs are shortened to <tmp>/ and the repo path to <COAL>/ so the
script stays readable on stage.
"""
from __future__ import annotations

import json
import os
import re
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from coal_app import core, tools  # noqa: E402


def tidy(obj):
    text = json.dumps(obj, indent=2) if not isinstance(obj, str) else obj
    text = text.replace(ROOT.replace("\\", "\\\\"), "<COAL>").replace(ROOT, "<COAL>")
    tmp = tempfile.gettempdir()
    for runs in (core.RUNS_DIR.replace("\\", "\\\\"), core.RUNS_DIR):
        text = re.sub(re.escape(runs) + r"[\\/]+run_[A-Za-z0-9_]+", "<tmp>", text)
    text = re.sub(re.escape(tmp.replace("\\", "\\\\")) + r"[\\/]+coal_[A-Za-z0-9_]+", "<tmp>", text)
    text = re.sub(re.escape(tmp) + r"[\\/]+coal_[A-Za-z0-9_]+", "<tmp>", text)
    return text.replace("\\\\", "/")


def call(fn, **kw):
    # The tool returns a compact summary itself, so nothing is dropped here.
    return fn(**kw)


def main():
    import shutil
    core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_demo_runs_")  # private, deleted below
    try:
        _main()
    finally:
        shutil.rmtree(core.RUNS_DIR, ignore_errors=True)


def _main():
    broken = core.example_paths("one_broken_gene")
    mito = core.example_paths("mito")
    r1 = call(tools.codon_align, protein_alignment_path=broken[0], dna_fasta_path=broken[1])
    e1 = tools.explain_report(r1["report_path"])
    r2 = call(tools.codon_align, protein_alignment_path=mito[0], dna_fasta_path=mito[1])
    code = next((h["suggested_codontable"] for h in r2["hints"] if h.get("suggested_codontable")), 2)
    r2["problems"] = r2["problems"][:3] + ([f"... {len(r2['problems']) - 3} more"] if len(r2["problems"]) > 3 else [])
    r3 = call(tools.codon_align, protein_alignment_path=mito[0], dna_fasta_path=mito[1], codontable=code)

    md = f"""# MCP demo: "which of my orthologs failed, and why?" (60 seconds)

A stage script for showing COAL as a tool an AI assistant can use. Everything below the script is the
real output of the tool functions, regenerated with `python demo/make_mcp_demo.py`.

## Setup (before going on stage)

1. `pip install -e ".[mcp]"` in the COAL repo, in its own virtual environment (not the one with `".[gui]"`).
   On Windows, `desktop\\setup_mcp.bat` does this and prints the two registrations below, filled in.
2. Add COAL to the assistant. Claude Code (use the full path of that environment's Python if it is not the
   `python` on your PATH):

   ```
   claude mcp add --transport stdio coal -- python -m coal_app --mcp
   ```

   Claude Desktop (`claude_desktop_config.json`). Use the full path of the Python that has COAL installed
   (`python -c "import sys; print(sys.executable)"`). Give the assistant full paths to the files in both cases:

   ```json
   {{
     "mcpServers": {{
       "coal": {{
         "command": "/full/path/to/python",
         "args": ["-m", "coal_app", "--mcp"]
       }}
     }}
   }}
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
{tidy({"protein_alignment_path": broken[0], "dna_fasta_path": broken[1]})}
```
Result:
```json
{tidy(r1)}
```

### 2. `explain_report` on that report

```
{tidy(e1)}
```

### 3. `codon_align` on mito (default genetic code 1)

```json
{tidy(r2)}
```

### 4. `codon_align` on mito again with `codontable={code}` (from the hint)

```json
{tidy(r3)}
```
"""
    out = os.path.join(ROOT, "demo", "mcp_demo.md")
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(md)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
