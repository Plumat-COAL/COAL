# Run notes (feed into SUMMARY.md)

## Tier A

- Normalisation used for grading: line endings (CRLF/CR -> LF) and line wrapping only (truth files are
  unwrapped; all four tools wrap FASTA at 60). Headers, order, case and every character compared exactly.
- Negative control: a truth file with one base changed is graded MISMATCH.
- 412/412 TRUTH for p2n14, p2n16, coal-py (-j1, default) and coal-cpp (-j1, default). No Tier A case produces
  an *output* difference from a documented divergence, so Tier A does not exercise COAL's deliberate differences.

### table10 warnings (v14 only)

p2n14 emits spurious "C does not correspond to TGA" warnings on table10 fam001, fam003, fam004 and fam005
(4, 16, 10 and 8 warnings; `expected_warnings` = 0). This is v14's documented table-10 defect: it never
implements the Euplotid code, so TGA stays a stop codon instead of Cys (p2n16 `tests/divergences.tsv`
`codontable_10_own_probe`; COAL README "Deliberate differences", code 10). The codon output is still correct
(TRUTH), so this is a warnings-only difference.

Why fam002 has no warning: its aligned CDS contains **no TGA codon at all** (0 TGA codons in-frame), so
there is nothing for the defect to fire on. In the other four the warning count equals the number of
TGA-as-Cys codons exactly (4, 16, 10, 8). p2n16, coal-py and coal-cpp emit 0 warnings on all five.

## Environment

- UAC is disabled on this machine (EnableLUA=0) and the benchmark session runs at High integrity, so
  `RAMMap64 -Et` runs without an elevation prompt. Each cold run verifies the standby list actually emptied
  (Memory\Standby Cache * Bytes before/after, recorded in `cold_method`); RAMMap exits 0 even when it cannot.
- Background at Tier B start: Defender (MsMpEng) ~0.5 core, Opera ~0.16 core, Discord ~0.03 core; WSL VM
  stopped. See processes_at_tierB_start.txt.

## Tier C design

- coal-py is timed once per batch (jobs=default): every Tier C file is < 8 Mnt, so automatic -j resolves to
  1 and `-j1` and default are the same code path; that row serves both headline columns.
- Runs: 1 warm-up + 5 timed for batch sizes 1/10/100; 1 warm-up + 3 timed for 1000 and the full set
  (user decision, to bound runtime).
- (tool, case) pairs that failed grading are not timed. To keep file sets identical across tools, a case is
  sampled only if every timed config passed grading on it; excluded counts are logged.
- Native batch mode (mode 2): N/A for COAL (no multi-alignment mode; several DNA files = one alignment's DNA
  split across files). Supplementary, non-headline "in-process API loop": one interpreter calling
  `main(argv)` per file, for coal-py and p2n16.
