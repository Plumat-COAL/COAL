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

## Tier B: p2n14 fails on the L35000 cases (genuine v14 limit)

- `n10_L35000` and `n100_L35000`: p2n14 exits 255 with empty stdout and the stderr
  `Too many capture groups (limit is 65535) in regex m/((U|T|C|Y|A)(U|T)G)(G(U|T).)... at pal2nal.pl line 1790.`
  This happened in all warm and cold runs.
- Cause: for each sequence v14 builds one regex by appending a reverse-translation pattern per residue (1 capture
  group for A/G/P/T, up to 7 for L/R/S; pal2nal.pl ~1759-1782), then runs `$nuc =~ /$qcodon/i` (line 1790).
  Perl caps capture groups per regex (the error says 65535; 65,533 is the most that compiles). These proteins need
  96,000-102,506 groups (every sequence in both cases is over the limit). The passing n1000_L10000 case peaks at
  28,254.
- What it means for v14: the limit depends on composition. For typical composition (about 2.8-2.9 groups per
  residue) it starts failing at roughly 22,000-24,000 aa. For a protein made only of L/R/S it could fail as low as
  about 9,400 aa. Tier A long_genes (up to 20,617 columns) all pass.
- Not a harness, Windows or Strawberry issue: reproduced outside the harness (plain `perl pal2nal.pl ... > out
  2> err`) under Strawberry Perl 5.42.3 (MSWin32), Git's Cygwin Perl 5.42.3, and Linux Perl 5.40.1 (WSL Ubuntu).
  All three give the same exit code (255), the same stderr, and the same 65,533 limit (probe script counts the
  groups). Timed runs: one ERROR row each, later runs SKIPPED. No median for p2n14 on these cases.

## Tier C: selenocysteine (documented v14 defect, not a harness issue)

- On SELD and FDHF (table11_large) p2n14 prints "ERROR: inconsistency between the following pep and nuc seqs",
  exits **0** and writes no output. v14 accepts `U` but no codon table maps it, so the pattern is one codon short.
  This is documented: grp-bork `tests/divergences.tsv` `fix_selenocysteine`, COAL README "U (selenocysteine)".
  The other tools' MISMATCH on those two cases only means "no reference": p2n14 gave nothing to compare
  against. These are to be reclassified KNOWN_DIV during triage.

## harness_commit history

- baseline rows: `b146167e713c`; warm Tier B: `ab92b8e18005`; cold Tier B and Tier C grading:
  `ab92b8e18005-dirty`. The "-dirty" is a false positive: the check also counted result files written under
  `benchmarks/results/`. `git diff ab92b8e -- benchmarks ':!benchmarks/results'` was empty at the time, so the
  code was identical. Fixed so result files are excluded from the check; rows keep the label as recorded.
