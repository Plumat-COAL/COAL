# COAL: likely judge questions, short honest answers

Judges: Apheris, PDBe, the Ellison Institute, CuspAI, Inherent, UCL. Rule: if a feature is marked
**[after merge]**, say "being merged now", or say nothing about it if the engine branch isn't in.

**1. How do you know it's correct on my data?**
We compared COAL with PAL2NAL v14's real output on 38 hand-made inputs covering every behaviour we found
(UTRs, frameshift codes, all 17 genetic codes, ambiguous bases, ID mismatches): 165 checks, 143
byte-identical, 12 rejected by both, 10 deliberate documented differences. On your data, COAL checks every
codon against its amino acid and tells you per gene if it doesn't fit, with the column, so problems can't
pass silently. Honest limit: the speed benchmarks use synthetic data, not real genomes, and the reference
inputs are small and hand-made.

**2. Why does 100x matter if PAL2NAL takes minutes?**
For one alignment it doesn't; small inputs are equally fast (about 0.05 s, measured). It matters at
tree-of-life scale: pipelines call this step once per gene family, thousands to millions of times, and
PAL2NAL slows down steeply on long genes (measured: 24.5 s vs 0.2 s for 100 genes x 10,000 codons). Half a
million sequences take COAL about 10 s. Speed is the second point; the first is that failures are explained
instead of silent.

**3. What happens with messy input?**
Built today: wrong genetic code -> COAL tests other codes and suggests the one that fits (one-click re-run
in the GUI). A frameshift -> the gene and column are named and shown in the alignment view. Invalid
characters -> a clear error (PAL2NAL silently drops a codon). Protein IDs like AT1G01010.1 vs DNA IDs
AT1G01010 -> COAL explains the mismatch. Errors give a non-zero exit code. **[after merge]** IDs matched
after removing version suffixes, and `--keep-going` so one bad gene doesn't stop the rest.

**4. Why not just use PAL2NAL?**
It works, and we treat it as the reference. But it always exits 0 even on errors, prints cryptic messages,
needs Perl on a Linux-style command line, and has no interface for bench scientists or AI agents. COAL
accepts PAL2NAL's options (there is even a drop-in `pal2nal` command), so a pipeline can switch without
being rewritten. PAL2NAL is GPL; COAL is an independent MIT implementation and contains none of its code.

**5. What about dN/dS?**
COAL produces the input that PAML (codeml), HyPhy and KaKs_Calculator need, and deliberately does not
estimate dN/dS itself: those tools are established and well validated, and re-implementing them would add
risk, not value. `-nogap` removes gap and stop-codon columns, as many dN/dS setups expect;
**[after merge]** PAML output format directly.

**6. Is this AI-generated code you can't maintain?**
Much of it was written with AI assistance (Claude), and we say so. What makes it maintainable is the
same as for any code: about 1,000 lines of dependency-free Python in the engine and 1,800 in the app
layer, a differential test suite against PAL2NAL's real output, 50+ app tests, CI configured for Windows,
macOS and Linux, and documented decisions (every deliberate difference from PAL2NAL is written down with its test).
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
CLUSTAL/PAML/codon output formats and `-blockonly` (output formats **[after merge]**); `-blockonly` and
`-html` give clear errors. The web GUI is new and has only been run on *[team: the laptops it was tested on]*; the Space is not deployed. The DIAMOND
preset has only been tested with a stand-in, not real genomes. The C++ core is **[after merge]**.
