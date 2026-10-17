# COAL: 3-minute pitch ("teach the room")

Rule for the speaker: say only what is on this page. Lines marked **[after merge]** describe engine features
being built on the engine branch (output formats, `--keep-going`, ID normalisation, C++ core). Say "this is
landing now", never "it does". If the engine branch is not merged by Sunday 10:00, skip those lines.

**Before going on:** laptop on mains, browser tabs open in this order: (1) `pitch/explainer.html`,
(2) the COAL GUI at http://127.0.0.1:7860 (`python -m coal_app`), (3) `pitch/backup/` folder. Zoom the
browser to 125-150 %. Run "One broken gene" once beforehand so nothing is cold.

## 0:00-0:30 Problem: what a codon alignment is for (tab 1, explainer plays itself, 20 s)

- "To ask how a gene changed through evolution, for example which genes were under selection in bats
  or in a pathogen, biologists line up the same gene from many species. Proteins are easy to line up.
  But the evidence for selection is in the DNA, three letters per amino acid."
- As the letters expand: "So we put the DNA back under the aligned protein, three letters per letter.
  Gaps become three gaps."
- As TGA flashes: "And when the DNA doesn't spell the protein, the result is wrong. Every downstream
  answer depends on catching this."

## 0:30-0:55 What wasn't working

- "The standard tool for this, PAL2NAL, is 20 years old and still used in pipelines everywhere."
- "Pipelines call it once per gene family, thousands of times, and its run time grows steeply with gene
  length: five times longer genes took 18 times longer in our measurements."
- "When something is wrong it prints a cryptic line, or nothing, and still reports success: its exit
  code is always 0. In one case we found it silently drops a codon."
- "It's a Perl script for the Linux command line. No interface for a bench scientist, and no way for an
  AI assistant to use it safely."

## 0:55-2:05 Live demo (tab 2, the COAL GUI)

| Time | Click | Say |
|---|---|---|
| 0:55 | (show the empty page) | "COAL runs locally, on this laptop. Your data stays here." |
| 1:00 | **Try an example: One broken gene** | "Twenty genes from twenty species. One has a single missing DNA letter." |
| 1:08 | Point at the red-orange row in the health table | "COAL names the gene and says in plain English what is wrong and where: a missing base at column 102. Every status has an icon and a word, so it never relies on colour alone." |
| 1:20 | Click **column 102** | "One click takes you to the exact codon. You can see the reading frame slip." |
| 1:30 | Scroll up, **Try an example: Mitochondrial genes** | "Now a classic mistake: mitochondrial genes run with the standard genetic code." |
| 1:38 | Point at the tip banner | "COAL notices that every gene fits the vertebrate mitochondrial code and offers to fix it." |
| 1:45 | Click **Re-run with genetic code 2** | "One click. Everything matches." |
| 1:52 | Point at the methods box, click copy | "And it writes the methods sentence for your paper, plus the exact command for your pipeline, and the PAL2NAL equivalent." |
| 2:00 | (one sentence) | "AI assistants get the same thing through MCP: two tools, `codon_align` and `explain_report`." |

On the current engine the broken-gene run stops at the bad gene ("COAL stopped: no alignment was written")
and the other 19 show as "not checked". **[after merge]** With `--keep-going` the same click writes the 19
good genes and flags the one bad one; say "the other 19 are aligned" only if you see it on screen.

## 2:05-2:35 Proof

- "Correctness: we test COAL against PAL2NAL's own output. 165 checks on 38 hand-made inputs and every
  option: 143 give byte-identical output, 12 are inputs both tools reject, and 10 are deliberate,
  documented differences where we think PAL2NAL is wrong, such as that dropped codon."
- "Speed, measured on one Windows laptop, end to end: 100 genes of 10,000 codons, PAL2NAL 24.5 seconds,
  COAL 0.2 seconds, 119 times faster. Ten thousand genes: 25 seconds down to under one second.
  Half a million genes, nearly a gigabyte: 10 seconds. For tiny inputs both take about 0.05 seconds; we
  don't claim a speed-up there."
- "The engine is pure Python with no dependencies, MIT licence. Continuous integration is set up for Windows, macOS and
  Linux." (Only say "passes on all three" if the GitHub Actions run is green when you present.)
- **[after merge]** "A C++ core for the largest jobs."

## 2:35-3:00 Vision

- "Tree-of-life biology, the kind of questions Hajk Drost asks, compares thousands of genes across
  hundreds of genomes. That means millions of codon alignments, and pipelines that call this step once per
  gene."
- "So we wrote a preset that goes from two species' raw files to codon alignments in one command: DIAMOND
  finds the matching genes, and COAL aligns all of them in a single run instead of once per gene." (Built and tested with a stand-in; say "a
  preset", not "we ran it on genomes".)
- "Old tools don't need replacing because they're wrong. They need to be fast, explain themselves, and
  be usable by a scientist, a pipeline, or an AI. That's COAL."

## Backup plan if the live demo fails

1. **Any GUI problem (won't start, browser error):** don't debug on stage. Switch to tab 3 and show
   `pitch/backup/` in order, keeping the same words and timings:
   - `1_broken_gene_frameshift.png`: the failed gene, its plain-English explanation and the frameshift column.
   - `2_mito_code1_hint.png`: the genetic-code tip.
   - `3_mito_tga_mismatch.png`: the TGA codons flagged.
   - `0_explainer_mismatch.png`: if the explainer won't play.
   These were rendered from the real viewer and tool output on the current engine (no GUI chrome).
2. **Laptop-only fallback, no Gradio needed:**
   `python -m coal_app.viewer examples/one_broken_gene/protein.aln.fasta --dna examples/one_broken_gene/dna.fasta -o view.html`,
   then open `view.html`.
3. **Team TODO before Sunday:** record a 70-second screen capture of the demo table above on the
   presenting laptop (the GUI could not be run in the overnight environment) and put it in `pitch/backup/`.
