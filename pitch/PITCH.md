# COAL: 3-minute pitch ("teach the room")

Rule for the speaker: say only what is on this page. `--keep-going`, safe ID pairing and the four output
formats are merged and tested in the Python engine (MRG-1..MRG-10, DEC-4), and the C++ engine (`coal-cpp`) mirrors
the merged code: 285/285 golden checks and 3,000/3,000 random C++ vs Python comparisons identical (overnight/core/LOG.md,
run 5 step 3). <!-- CHANGED(PITCH-2): was "Lines marked **[C++ core]** depend on the C++ engine passing its comparison tests against the merged code (Person A's step 3): if it hasn't by Sunday 10:00, skip them." Step 3 passed at 47d2168 --> <!-- CHANGED(PITCH-1): was "[after the C++/formats branch is merged]" for --keep-going, ID normalisation and the C++ core; the Python features are merged now --> <!-- CHANGED(DEC-4): was "output formats" on the teammate's branch and "`coal --help` still lists only -output fasta" -->
<!-- CHANGED(R2-D12): was "[after merge]", which read as the (already merged) engine-fixes branch -->

**Before going on:** laptop on mains, browser tabs open in this order: (1) `pitch/explainer.html`,
(2) the COAL GUI at http://127.0.0.1:7860 (`python -m coal_app`), (3) `pitch/backup/` folder. Zoom the
browser to 125-150 %. Run "One broken gene" once beforehand so nothing is cold.

## 0:00-0:30 Problem: what a codon alignment is for (tab 1, explainer plays itself: last step at about 16 s, then it holds)

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
- "It's a Perl script, or a web form that handles one alignment at a time and can't be put in a pipeline.
  Neither explains what went wrong to a bench scientist, and there's no way for an AI assistant to use it safely."

## 0:55-2:05 Live demo (tab 2, the COAL GUI)

| Time | Click | Say |
|---|---|---|
| 0:55 | (show the empty page) | "COAL runs locally, on this laptop. Your data stays here." |
| 1:00 | **Try an example: One broken gene** | "Twenty genes from twenty species. One has a single missing DNA letter." |
| 1:08 | Point at the red-orange row in the health table | "COAL names the gene and says in plain English what is wrong, and the app estimates where: probably a missing base at column 102. Every status has an icon and a word, so it never relies on colour alone." |
| 1:20 | Click **column 102** | "One click takes you to the exact codon where the frame shifts. The viewer re-reads the DNA from there in the new frame, so the codons after it line up again." |
| 1:30 | Scroll up, **Try an example: Mitochondrial genes** | "Now a classic mistake: mitochondrial genes run with the standard genetic code." |
| 1:38 | Point at the tip banner | "COAL notices that every gene fits the vertebrate mitochondrial code and offers to fix it." |
| 1:45 | Click **Re-run with genetic code 2** | "One click. Everything matches." |
| 1:52 | Point at the methods box, click copy | "And it writes the methods sentence for your paper, plus the exact command for your pipeline, and the closest PAL2NAL command." <!-- CHANGED(FIN-A5) --> |
| 2:00 | (one sentence) | "AI assistants get the same thing through MCP: two tools, `codon_align` and `explain_report`." |

The app runs with `--keep-going`, so this click writes the 19 good genes and flags the one bad one (checked on the
merged code: 19 OK, 1 failed, alignment written). You can say "the other 19 are aligned"; check it on screen during the
rehearsal. <!-- CHANGED(PITCH-1): was "stops at the bad gene" on the pre-merge engine, with --keep-going as a future feature -->

## 2:05-2:35 Proof

<!-- CHANGED(R2-D2): was 143 / 94 / 10 (stale after the engine merge); now the harness's current counts, guarded by tests_app/test_docs_numbers.py. See overnight/integration/CHANGES.md -->
- "Correctness: we test COAL against PAL2NAL version 14's own output: 285 checks on 39 hand-made inputs, in
  all four output formats. 238 compare the output byte for byte, and 234 of those also compare the warnings.
  24 check that COAL rejects the inputs PAL2NAL fails on. 23 are documented, deliberate differences, such as that dropped codon, newer genetic-code tables or `-nogap` using
  the chosen genetic code; for those we check COAL's own behaviour. <!-- CHANGED(FIN-C3): was "where we think PAL2NAL is wrong" (not true of every one) --> And on every
  change, CI runs both tools on 3,000 random inputs in every format. In short: COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate differences (newer NCBI genetic-code tables for codes 3, 10, 11 and 13; `-nogap` removing the selected code's stop codons; the selected code's start codons for the first Met; pairing by IDs that match only after normalising where PAL2NAL v14's file order would give a protein another gene's DNA; and inputs where PAL2NAL v14 writes malformed output or fails). Tested: all 285 stored reference checks pass (238 byte-identical, 4 of them with deliberately different warnings; 24 where both tools fail; 23 deliberate output or exit-status differences: 6 deliberate errors, 3 different outputs, 14 `-nogap` checks with a non-standard code), and in 20,000 random comparisons with the real PAL2NAL v14 (5,000 cases x 4 output formats) 13,826 were byte-identical, 4,681 both failed and 1,493 were deliberate differences, each checked to be exactly the documented change; 0 unexplained differences." <!-- CHANGED(FIN-C2/FIN-C3/FIN-C10/FIN-C20): was the "byte-identical ... (well-formed input)" claim, "every difference is one where PAL2NAL silently produces wrong output" and a 6,000-case run not reproducible from the repo; now the documented exceptions and the FIN-C10 fuzz numbers --> <!-- CHANGED(FZ-claim) --> <!-- CHANGED(DEC-4): was 165 / 147 / 141 / 12 / 6, FASTA only; fuzz sentence new --> <!-- CHANGED(FZ-5): 246/240/9 -> 232/228/23 -->
  **Spoken version (say this; the paragraph above is the reference if a judge asks for detail):** "285 checks
  against PAL2NAL's own output, all pass. 20,000 random comparisons with the real PAL2NAL: zero unexplained
  differences. Where we differ on purpose, it's written down and tested." <!-- CHANGED(PITCH-2): new; 285 and 20,000 as in the paragraph above (tests_app/test_docs_numbers.py, README Compatibility) -->
- "Speed, measured on one Windows laptop, end to end: 100 genes of 10,000 codons, PAL2NAL 24.5 seconds,
  COAL 0.2 seconds, 119 times faster. Ten thousand genes: 25 seconds down to under one second.
  Half a million genes, nearly a gigabyte: 10 seconds with 4 processes, 25 on one. For tiny inputs both take about 0.05 seconds; we
  don't claim a speed-up there."
  (Python engine, benchmarks/results/bench_windows_20261016_2356: measured before the weekend's fixes and the merge.
  Person A's step-4 re-run on the merged code (benchmarks/results/cpp_vs_python_linux_4084d1f.md) compared only
  Python with C++, not with PAL2NAL. A pre-submission check against the real PAL2NAL v14 on a 2-core Linux VM gave
  100 x 10,000 codons: PAL2NAL 41.4 s, COAL Python 0.28 s, C++ 0.038 s; 10,000 x 500: PAL2NAL 31.3 s, Python
  0.99-1.6 s, C++ 0.21 s. These are **provisional and not in the repo: re-measure on the presenting laptop** before
  quoting any of them. <!-- PITCH-1 --> <!-- CHANGED(PITCH-2): was "Swap in the numbers from Person A's step-4 re-run on the merged code if they differ." -->)
- "The engine is pure Python with no dependencies, MIT licence. Continuous integration is set up for Windows, macOS and
  Linux." (Only say "passes on all three" if the GitHub Actions run is green when you present.)
- "There's also a C++ engine for the largest jobs. It passes all 285 checks, and in 3,000 random comparisons its
  output is identical to our Python engine's. It's three to six times faster than the Python engine on the
  same machine." (benchmarks/results/cpp_vs_python_linux_4084d1f.md, 2-core Linux VM: 5.3-6.3x on 100 x 10,000 and 10,000 x 500,
  3.1-4.0x on 100,000 x 500. It is built from source with CMake; the pip package is the Python engine only.) <!-- CHANGED(PITCH-2): was "**[C++ core]** "A C++ core for the largest jobs." (Speed numbers above are for the Python engine; quote C++ numbers only from Person A's step-4 benchmark on the merged code.)" -->

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
