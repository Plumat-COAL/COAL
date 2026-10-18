# COAL: 3-minute pitch ("teach the room")

Rule for the speaker: say only what is on this page. `--keep-going`, safe ID pairing and the four output
formats are merged and tested in the Python engine (MRG-1..MRG-10, DEC-4). The C++ engine is merged too
(it matches the Python engine on 285/285 stored checks and on thousands of random inputs). <!-- CHANGED(PITCH-1): was "[after the C++/formats branch is merged]" for --keep-going, ID normalisation and the C++ core; the Python features are merged now --> <!-- CHANGED(DEC-4): was "output formats" on the teammate's branch and "`coal --help` still lists only -output fasta" -->
<!-- CHANGED(R2-D12): was "[after merge]", which read as the (already merged) engine-fixes branch -->

**Before going on:** laptop on mains, browser tabs open in this order: (1) `site/explainer.html`,
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
<!-- CHANGED(REV-2): was 120 spoken words (about 48 s) for a 25 s slot -->

- "The standard tool for this, PAL2NAL, is 20 years old and still in pipelines everywhere. Its run time grows
  steeply with gene length: five times longer genes took 18 times longer in our measurements."
- "And when something is wrong it prints a cryptic line, or nothing, and still exits with success. In one case we
  found it silently drops a codon."

## 0:55-2:05 Live demo (tab 2, the COAL GUI)

| Time | Click | Say |
|---|---|---|
| 0:55 | (show the empty page) | "COAL runs locally, on this laptop. Your data stays here." |
| 1:00 | **Try an example: One broken gene** | "Twenty synthetic genes, one per species. One has a single missing DNA letter." <!-- CHANGED(REV-2): the example data are synthetic --> |
| 1:08 | Point at the red row in the health table | "COAL names the gene and says in plain English what is wrong, and the app estimates where: probably a missing base at column 102. Every status has an icon and a word, so it never relies on colour alone." |
| 1:20 | Click **column 102** | "One click takes you to the exact codon where the frame shifts. The viewer re-reads the DNA from there in the new frame, so the codons after it line up again." |
| 1:30 | Scroll up, **Try an example: Mitochondrial (wrong code)** <!-- CHANGED(AMB-3): was "Mitochondrial genes" --> | "Now a classic mistake: mitochondrial genes run with the standard genetic code." |
| 1:38 | Point at the banner "Wrong genetic code?" <!-- CHANGED(AMB-4) --> | "COAL notices that every gene fits the vertebrate mitochondrial code and offers to fix it." |
| 1:45 | Click **Re-run with genetic code 2** | "One click. Everything matches." |
| 1:52 | Point at the methods box, click copy | "And it writes the methods sentence for your paper, plus the exact command for your pipeline, and the closest PAL2NAL command." <!-- CHANGED(FIN-A5) --> |
| 2:00 | (one sentence) | "AI assistants get the same thing through MCP: two tools, `codon_align` and `explain_report`." |

The app runs with `--keep-going`, so this click writes the 19 good genes and flags the one bad one (checked on the
merged code: 19 OK, 1 failed, alignment written). You can say "the other 19 are aligned"; check it on screen during the
rehearsal. <!-- CHANGED(PITCH-1): was "stops at the bad gene" on the pre-merge engine, with --keep-going as a future feature -->

## 2:05-2:35 Proof
<!-- CHANGED(REV-2): the spoken part was 334 words (over two minutes) because the README's full compatibility claim
was quoted; now about 110 spoken words, with the full claim kept below as reference for questions -->

- "Correctness: we test against PAL2NAL's own output. 285 stored checks in all four formats, plus 20,000 random
  comparisons with the real PAL2NAL v14: zero unexplained differences. Every difference is documented and
  deliberate, like that dropped codon."
- "We built this in a weekend with Claude. The reason you can trust it is that the real PAL2NAL checked every
  change." <!-- REV-2: how AI was used; QA.md Q6 has the detail -->
- "Speed: 100 long genes, PAL2NAL 24 seconds, COAL 0.2: over a hundred times faster. And a C++ core, which we
  measured at three to six times faster again." (Python numbers: benchmarks/results/bench_windows_20261016_2356,
  Windows, before the weekend's fixes; C++ vs Python: benchmarks/results/cpp_vs_python_linux_4084d1f.md, merged
  code, 2-core Linux VM, 3.1-6.3x. Re-run run_bench.py on the presenting laptop if time allows.)
- "A Python engine with no required dependencies, an optional C++ engine for the biggest jobs, MIT licence." Add "CI passes on Windows, macOS and Linux" only if the GitHub
  Actions run is green when you present.

Reference, not read aloud (for questions; guarded by tests_app/test_docs_numbers.py):


<!-- CHANGED(R2-D2): was 143 / 94 / 10 (stale after the engine merge); now the harness's current counts, guarded by tests_app/test_docs_numbers.py. See overnight/integration/CHANGES.md -->
- Correctness: we test COAL against PAL2NAL version 14's own output: 285 checks on 39 hand-made inputs, in
  all four output formats. 238 compare the output byte for byte, and 234 of those also compare the warnings.
  24 check that COAL rejects the inputs PAL2NAL fails on. 23 are documented, deliberate differences, such as that dropped codon, newer genetic-code tables or `-nogap` using
  the chosen genetic code; for those we check COAL's own behaviour. <!-- CHANGED(FIN-C3): was "where we think PAL2NAL is wrong" (not true of every one) --> And on every
  change, CI runs both tools on 3,000 random inputs in every format. In short: COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate differences (newer NCBI genetic-code tables for codes 3, 10, 11 and 13; `-nogap` removing the selected code's stop codons; the selected code's start codons for the first Met; pairing by IDs that match only after normalising where PAL2NAL v14's file order would give a protein another gene's DNA; and inputs where PAL2NAL v14 writes malformed output or fails). Tested: all 285 stored reference checks pass (238 byte-identical, 4 of them with deliberately different warnings; 24 where both tools fail; 23 deliberate output or exit-status differences: 6 deliberate errors, 3 different outputs, 14 `-nogap` checks with a non-standard code), and in 20,000 random comparisons with the real PAL2NAL v14 (5,000 cases x 4 output formats) 13,826 were byte-identical, 4,681 both failed and 1,493 were deliberate differences, each checked to be exactly the documented change; 0 unexplained differences." <!-- CHANGED(FIN-C2/FIN-C3/FIN-C10/FIN-C20): was the "byte-identical ... (well-formed input)" claim, "every difference is one where PAL2NAL silently produces wrong output" and a 6,000-case run not reproducible from the repo; now the documented exceptions and the FIN-C10 fuzz numbers --> <!-- CHANGED(FZ-claim) --> <!-- CHANGED(DEC-4): was 165 / 147 / 141 / 12 / 6, FASTA only; fuzz sentence new --> <!-- CHANGED(FZ-5): 246/240/9 -> 232/228/23 -->
- Speed detail: 100 genes of 10,000 codons, PAL2NAL 24.5 s, COAL 0.2 s (119x). Ten thousand genes: 25 s down to
  under one second. Half a million genes, nearly a gigabyte: 10 s with 4 processes, 25 s on one. For tiny inputs
  both take about 0.05 s; we don't claim a speed-up there.

## 2:35-3:00 Vision
<!-- CHANGED(REV-2): was 105 spoken words for a 25 s slot -->

- "Comparative genomics now aligns thousands of genes across hundreds of genomes, so we added a preset that goes from
  two species' raw files to codon alignments in one command." (Built and tested with a stand-in: say "a preset",
  not "we ran it on genomes".)
- "Old tools don't need replacing because they're wrong. They need to be fast, explain themselves, and be usable by
  a scientist, a pipeline, or an AI. That's COAL."

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
