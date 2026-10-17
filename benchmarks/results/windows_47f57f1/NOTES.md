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

## GitHub issue draft: coal-py commits ~630 MB per process once it imports numpy (default -j: up to ~11 GB)

**Title:** coal-py: each worker commits ~630 MB after the lazy numpy import, so default `-j` memory scales with core
count, not input size

**Summary.** On Windows (i5-13600KF, 20 threads, numpy 2.4.6 + bundled `libscipy_openblas64_`), coal-py's
peak committed memory jumps by about 630 MB in every process that imports numpy. At the default `-j` on 100k
sequences that means 9.0-10.8 GB, against 0.69 GB at `-j1`, for a 1.9x speedup. None of it is input data.

**Measurements** (Tier B, `runs.csv`, peak committed memory of the whole Job Object):

| case | codons | -j1 | default -j (workers) | -j12 |
|---|---|---|---|---|
| n1000_L400 | 0.4 M | 18.6 MB | 18.7 MB (1) | 18.6 MB (1) |
| n10000_L400 | 4.0 M | 647 MB | 98 MB (3) | 99 MB (3) |
| n100000_L400 | 39.4 M | 688 MB | 9,024-10,849 MB (20) | 7,746 MB (12) |

**Where it comes from** (code reading only, COAL `47f57f1`):

- `core.py:353-360` `_quick_offset0`: for sequences of at least `NP_MIN_CODONS` = 300 codons, a per-process
  counter accumulates codons. Once it reaches `NP_AFTER_CODONS` = 2,000,000 (`core.py:371`), `_numpy()` imports numpy
  (`core.py:374-387`). numpy is used only for a uint8 table lookup and compare (`_np_quick`, `core.py:396`), never
  for BLAS.
- So the import happens in each process that individually checks 2 M codons, and every number above matches:
  - n1k: no process reaches 2 M.
  - n10k at -j1: one process, 4 M codons, imports.
  - n10k at default (3 workers): about 1.3 M each, none imports.
  - n100k at -j12: about 3.3 M each, all 12 import (7,746 MB is about 12 × 636 MB).
  - n100k at default (20 workers): about 1.97 M each, right at the threshold, so some workers import and some don't.
    That is why that config's peak varies from run to run (9.0-10.8 GB).
- Things ruled out: work batches are about 2 Mnt (`BATCH_NT`, `pipeline.py:52`) with at most 2×jobs in flight
  (`pipeline.py:307-323`). `DnaIndex` keeps only offsets and reads sequences on demand (`seqio.py:224-303`). Output is
  streamed through a temp file (`pipeline.py:638`, `693-698`). None of these accounts for ~630 MB per process.
- Likely mechanism (**hypothesis, not yet verified**): numpy's bundled OpenBLAS allocates per-thread buffers when it
  loads, one thread per logical CPU by default (20 here, about 20 × 32 MB). On Windows this is committed memory.
  If so, the per-process cost grows with core count, and so does the number of workers, so total memory would grow
  roughly with cores² on bigger machines.
- To verify (after the timed tiers, untimed): peak commit of `python -c "import numpy"` with and without
  `OPENBLAS_NUM_THREADS=1`. Also coal-py n100k -j12 with `OPENBLAS_NUM_THREADS=1`, and again with `COAL_NO_NUMPY=1`.

**Suggested fixes, in order:**

1. **Limit BLAS threads before the import.** In `_numpy()`, `os.environ.setdefault(...)` for
   `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `MKL_NUM_THREADS` = "1" before `import numpy`. COAL never calls
   BLAS, so this costs nothing. It must run before numpy is first imported in that process. A worker initializer
   (`_init_worker`) is the other place to do it.
2. **Make the threshold deterministic, or drop numpy here.** The per-process 2 M-codon trigger makes memory depend
   on how batches land on workers. An alternative is a numpy-free fast path (for example `bytes.translate` over a
   codon-index encoding), which would remove the dependency's start-up and memory cost entirely.
3. **Cap workers by memory as well as by size.** `effective_jobs` (`pipeline.py:100`) could also limit workers by
   available commit divided by a per-worker estimate. This only mitigates; 1 removes the cause.
- Chunked work or shared memory would **not** help: the per-worker memory is not input data (see "Things ruled out"
  above).

**Not changed in this benchmark run.** COAL is frozen at `47f57f1`.
