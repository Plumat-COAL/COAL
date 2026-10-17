# COAL: likely judge questions, short honest answers

Judges: Apheris, PDBe, the Ellison Institute, CuspAI, Inherent, UCL. Everything below is merged and tested, including
the C++ engine. <!-- CHANGED(REV-3): was a rule for features "[after the C++/formats branch is merged]"; all are merged now -->

**1. How do you know it's correct on my data?**
We compared COAL with PAL2NAL v14's real output on 39 hand-made inputs covering every behaviour we found
(UTRs, frameshift codes, all 17 genetic codes, ambiguous bases, ID mismatches): 285 checks on 39 hand-made
inputs, in all four output formats. 238 compare COAL's output with PAL2NAL's byte for byte, and 234 of those also
compare the warnings; 24 check that COAL rejects the inputs PAL2NAL fails on (exit code only, not the message); 23 are
documented, deliberate differences (14 of them `-nogap` with another genetic code: COAL removes that code's stop codons) <!-- CHANGED(FZ-5): was 246/240/9 -->, where we check COAL's own behaviour instead. On every change, CI also runs both
tools on 3,000 random inputs x 4 formats (wrong DNA, masked genes, UTRs, frameshifts, every genetic code, -nogap,
-nomismatch, -blockonly): identical apart from the documented differences. <!-- CHANGED(DEC-4): was 165/147/141/12/6 -->
In one sentence: COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate differences (newer NCBI genetic-code tables for codes 3, 10, 11 and 13; `-nogap` removing the selected code's stop codons; the selected code's start codons for the first Met; pairing by IDs that match only after normalising where PAL2NAL v14's file order would give a protein another gene's DNA; and inputs where PAL2NAL v14 writes malformed output or fails). Tested: all 285 stored reference checks pass (238 byte-identical, 4 of them with deliberately different warnings; 24 where both tools fail; 23 deliberate output or exit-status differences: 6 deliberate errors, 3 different outputs, 14 `-nogap` checks with a non-standard code), and in 20,000 random comparisons with the real PAL2NAL v14 (5,000 cases x 4 output formats) 13,826 were byte-identical, 4,681 both failed and 1,493 were deliberate differences, each checked to be exactly the documented change; 0 unexplained differences. <!-- CHANGED(FIN-C2/FIN-C3/FIN-C10/FIN-C20): was the "byte-identical ... (well-formed input)" claim, "every difference is one where PAL2NAL silently produces wrong output" and a 6,000-case run not reproducible from the repo; now the documented exceptions and the FIN-C10 fuzz numbers --> <!-- CHANGED(FZ-claim) -->
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
output as PAL2NAL v14 where PAL2NAL writes one (it refuses only some short or frameshifted genes, with
"inconsistency"), plus <!-- CHANGED(FIN-C17): was "(it never refuses)" --> a plain warning naming the gene, and the GUI marks it failed;
`--strict` makes it an error instead. <!-- CHANGED(DEC-1): new sentence --> Protein IDs like AT1G01010.1 vs DNA IDs
AT1G01010 -> COAL explains the mismatch. Errors give a non-zero exit code. IDs that match only
after removing version suffixes are paired by name, with a notice,, and `--keep-going` so one bad gene doesn't stop the rest.

**4. Why not just use PAL2NAL?**
It works, and we treat it as the reference. But it always exits 0 even on errors, prints cryptic messages,
is a Perl script (or a web form that takes one alignment at a time and can't be used from a pipeline), and
explains nothing to bench scientists or AI agents. COAL
accepts PAL2NAL's options and writes all four of its output formats as PAL2NAL v14 does, byte for byte apart from
the documented differences, and the drop-in `pal2nal` <!-- CHANGED(FIN-C3): was "byte for byte" without the exceptions -->
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
same as for any code: about 3,100 lines of dependency-free Python in the engine and 4,600 in the app
layer, a differential test suite against PAL2NAL's real output, more than 300 app tests, CI on Windows, macOS
and Linux, including the GUI tests against real Gradio and the C++ engine on all three (open the latest GitHub
Actions run on main before you say "green") <!-- CHANGED(SHIP-5): was "all green in run 39006043852"; later runs were red on Windows until SHIP-2 -->, and documented decisions (every deliberate difference from PAL2NAL is written down with its test).
<!-- CHANGED(FINAL-2): was "the Windows job was red until the round-2 fixes" (it was still red after round 2); now cites green run 39006043852. -->
<!-- CHANGED(R2-D2, R2-D3): was "about 1,050 lines ... 2,150 ... 115 app tests ... green on the app branch; the new GUI and MCP jobs have not run yet" (both false). See overnight/integration/CHANGES.md -->
<!-- CHANGED(R3-D4b): app test count 140 -> 170 (round 3 added tests; the docs guard requires it). See overnight/integration/CHANGES.md -->
<!-- CHANGED(DEC-4): engine line count 1,300 -> 1,550 (output formats; the docs guard requires it). See overnight/integration/CHANGES.md -->
<!-- CHANGED(MRG-2..4): engine line count 1,550 -> 3,100 (run report, keep-going, resume, progress from the merged branch) -->
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
output formats and -blockonly" --> The desktop app has been run on Windows 11 (the command line byte-identical to Linux) and in CI on Linux; the Space is not deployed. The DIAMOND
preset has only been tested with a stand-in, not real genomes. The C++ core is merged and matches the Python engine (285/285 stored checks, thousands of random inputs); it is
built from source with CMake, not yet shipped in the installer.

**11. Why not just fix PAL2NAL?** <!-- REV-3: new -->
It's GPL Perl, and we wanted an MIT engine anyone can embed, plus a C++ core. We kept PAL2NAL as the reference
that judged us: every change was compared with its output.

**12. What if the gene really has a frameshift, e.g. a pseudogene?** <!-- REV-3: new -->
Then COAL is right to stop on it: a frameshifted gene has no valid codon alignment. With `--keep-going` it is left
out and named, so you can drop it or fix its sequence; the rest of the family is written.
