# COAL: likely judge questions, short honest answers

Judges: Apheris, PDBe, the Ellison Institute, CuspAI, Inherent, UCL. Rule: if a feature is marked
**[after the C++/formats branch is merged]**, say "being merged now", or say nothing about it if that branch
isn't in (the engine-fixes branch is already merged). <!-- CHANGED(R2-D12): was "[after merge]" -->

**1. How do you know it's correct on my data?**
We compared COAL with PAL2NAL v14's real output on 38 hand-made inputs covering every behaviour we found
(UTRs, frameshift codes, all 17 genetic codes, ambiguous bases, ID mismatches): 285 checks on 39 hand-made
inputs, in all four output formats. 238 compare COAL's output with PAL2NAL's byte for byte, and 234 of those also
compare the warnings; 24 check that COAL rejects the inputs PAL2NAL fails on (exit code only, not the message); 23 are
documented, deliberate differences (14 of them `-nogap` with another genetic code: COAL removes that code's stop codons) <!-- CHANGED(FZ-5): was 246/240/9 -->, where we check COAL's own behaviour instead. On every change, CI also runs both
tools on 3,000 random inputs x 4 formats (wrong DNA, masked genes, UTRs, frameshifts, every genetic code, -nogap,
-nomismatch, -blockonly): identical apart from the documented differences. <!-- CHANGED(DEC-4): was 165/147/141/12/6 -->
In one sentence: COAL is byte-identical to PAL2NAL v14 on well-formed input: all 285 stored reference checks pass (238 byte-identical, 24 where both tools fail, 23 documented deliberate differences), and in 20,000 random comparisons (5,000 cases x 4 output formats) and 6,000 independent stratified cases every difference from PAL2NAL is a documented, deliberate one, where PAL2NAL silently produces wrong or malformed output; 0 undocumented differences. <!-- CHANGED(FZ-claim) -->
On your data, COAL checks every
<!-- CHANGED(R2-D2): was 143 / 94 / 10; guarded by tests_app/test_docs_numbers.py -->
codon against its amino acid and tells you per gene if it doesn't fit, with the column, so problems can't
pass silently (with `-nomismatch` the command line, like PAL2NAL, prints no mismatch warnings; the GUI and
the AI tools re-check and still name the genes and count the removed columns).
<!-- CHANGED(R2-D5) --> Honest limit: the speed benchmarks use synthetic data, not real genomes, and the reference
inputs are small and hand-made.

**2. Why does 100x matter if PAL2NAL takes minutes?**
For one alignment it doesn't; small inputs are equally fast (about 0.05 s, measured). It matters at
tree-of-life scale: pipelines call this step once per gene family, thousands to millions of times, and
PAL2NAL slows down steeply on long genes (measured: 24.5 s vs 0.2 s for 100 genes x 10,000 codons). Half a
million sequences take COAL about 10 s with 4 processes (25 s on one). Speed is the second point; the first is that failures are explained
instead of silent.

**3. What happens with messy input?**
Built today: wrong genetic code -> COAL tests other codes and suggests the one that fits (one-click re-run
in the GUI). A frameshift -> the gene and column are named and shown in the alignment view. Invalid
characters -> a clear error (PAL2NAL silently drops a codon). DNA that doesn't belong to the protein -> the same
output as PAL2NAL (it never refuses), plus a plain warning naming the gene, and the GUI marks it failed;
`--strict` makes it an error instead. <!-- CHANGED(DEC-1): new sentence --> Protein IDs like AT1G01010.1 vs DNA IDs
AT1G01010 -> COAL explains the mismatch. Errors give a non-zero exit code. **[after the C++/formats branch is merged]** IDs matched
after removing version suffixes, and `--keep-going` so one bad gene doesn't stop the rest.

**4. Why not just use PAL2NAL?**
It works, and we treat it as the reference. But it always exits 0 even on errors, prints cryptic messages,
is a Perl script (or a web form that takes one alignment at a time and can't be used from a pipeline), and
explains nothing to bench scientists or AI agents. COAL
accepts PAL2NAL's options and writes all four of its output formats byte for byte, and the drop-in `pal2nal`
command defaults to clustal like PAL2NAL, so an unchanged PAL2NAL command line works: a pipeline switches by changing
the command name. <!-- CHANGED(DEC-4): was "pipelines that rely on PAL2NAL's default clustal output must add
-output fasta" --> PAL2NAL is GPL; COAL is an independent MIT implementation and contains none of its code.

**5. What about dN/dS?**
COAL produces the input that PAML (codeml), HyPhy and KaKs_Calculator need, and deliberately does not
estimate dN/dS itself: those tools are established and well validated, and re-implementing them would add
risk, not value. `-nogap` removes gap and stop-codon columns, as many dN/dS setups expect, and `-output paml`
writes codeml's input directly. <!-- CHANGED(DEC-4): PAML output was "[after the C++/formats branch is merged]" -->

**6. Is this AI-generated code you can't maintain?**
Much of it was written with AI assistance (Claude), and we say so. What makes it maintainable is the
same as for any code: about 3,000 lines of dependency-free Python in the engine and 2,600 in the app
layer, a differential test suite against PAL2NAL's real output, more than 200 app tests, CI on Windows, macOS
and Linux (all green in run 39006043852, including the GUI tests against real Gradio on Windows; still open the
latest GitHub Actions run before you say "green"), and documented decisions (every deliberate difference from PAL2NAL is written down with its test).
<!-- CHANGED(FINAL-2): was "the Windows job was red until the round-2 fixes" (it was still red after round 2); now cites green run 39006043852. -->
<!-- CHANGED(R2-D2, R2-D3): was "about 1,050 lines ... 2,150 ... 115 app tests ... green on the app branch; the new GUI and MCP jobs have not run yet" (both false). See overnight/integration/CHANGES.md -->
<!-- CHANGED(R3-D4b): app test count 140 -> 170 (round 3 added tests; the docs guard requires it). See overnight/integration/CHANGES.md -->
<!-- CHANGED(DEC-4): engine line count 1,300 -> 1,550 (output formats; the docs guard requires it). See overnight/integration/CHANGES.md -->
<!-- CHANGED(MRG-2..4): engine line count 1,550 -> 3,000 (run report, keep-going, resume, progress from the merged branch) -->
*[Team: say who reviewed which parts. Only claim review you actually did.]*

**7. How will this be maintained after the weekend?**
It's public on GitHub under MIT, with tests that run on every change, so others can safely contribute.
*[Team: replace this sentence with who will maintain it and for how long, e.g. "X will keep it for her
PhD work". Don't promise more than you will do.]*

**8. Does my data leave my machine? (Apheris)**
No. The engine, GUI and MCP server run locally; the GUI is set up to make no external calls (system fonts
only, Gradio analytics off). A Hugging Face Space is prepared for people who want a hosted demo, but it is not deployed and is
opt-in. Honest limit: we haven't audited Gradio's own network behaviour beyond turning analytics off.

**9. How do AI agents use it, and how do you stop them making things up?**
Through MCP: `codon_align` returns the output file, counts, hints and per-gene explanations as structured
data, and `explain_report` turns the JSON report into plain English. The agent reports what COAL found,
with the same columns and gene names a human sees, instead of guessing from raw error text. It can't
"fix" a sequence; it can only re-run with options COAL suggests, such as another genetic code.

**10. What doesn't it do yet?**
PAL2NAL's `-html` web-server mode (a clear error) and Gblocks-format input. <!-- CHANGED(DEC-4): was "CLUSTAL/PAML/codon
output formats and -blockonly" --> The web GUI is new and has only been run on *[team: the laptops it was tested on]*; the Space is not deployed. The DIAMOND
preset has only been tested with a stand-in, not real genomes. The C++ core is **[after the C++/formats branch is merged]**.
