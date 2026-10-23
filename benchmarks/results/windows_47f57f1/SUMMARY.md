# COAL Windows benchmark: summary

Tools: COAL `47f57f1` (coal-py, coal-cpp), PAL2NAL v14 (Perl, "p2n14"), PAL2NAL v16 grp-bork Python port
`da572f7` ("p2n16"). Machine: i5-13600KF (6 P-cores with SMT + 8 E-cores = 14 cores / 20 threads), 32 GB, NVMe,
Windows 11. Details: `env.json`. Raw data: `runs.csv`, `correctness.csv.gz` (merged; the per-tier `correctness_*.csv` and `diffs/` are on the `bench-windows-47f57f1` tag). Run log and design notes: `status.log`,
`NOTES.md`. **This run was cut short by a deadline; see "Scope reductions due to time" at the end.**

## Headlines

**The headline is single-threaded speed.** coal-cpp at `-j1` is 74–1,115× faster than p2n14 on the Tier B workloads.
Extra threads add little: **coal-cpp plateaus at about 1.2× from `-j2`** (Tier E). At 1M sequences, **p2n14
took 3,714 s (62 min) and 18.6 GB, output correct (TRUTH)**. That is a single run with no warm-up, made overnight.
So coal-cpp `-j1` (12.52 s) is **297× faster** (347× at default `-j`). **p2n16 did not finish (DNF, >900 s)**,
so against p2n16 coal-cpp `-j1` is **more than 71.9× faster**, a lower bound.

Speed-ups against p2n14, from warm medians:

| workload | coal-cpp -j1 | coal-cpp default (-j20) | coal-py -j1 | coal-py default |
|---|---|---|---|---|
| 100k seqs × 400 codons (p2n14 239 s) | **223×** | 270× | 47× | 90× |
| 1k seqs × 10,000 codons (p2n14 280 s) | **1,115×** | 1,296× | 292× | 325× |
| 1k seqs × 400 codons (p2n14 2.3 s) | 74× | 86× | 12.5× | 12.8× |
| 100 real NCBI families, one process per file (p2n14 125 ms/file) | **6.9×** (18 ms/file) | 6.9× | **0.68×** (184 ms/file) | 0.68× |
| 1M seqs × 400 codons, vs p2n16 DNF (>900 s) | >71.9× (12.5 s) | >84× (10.7 s) | >17.9× (50.3 s) | >34× (26.4 s, single run) |
| 1M seqs × 400 codons (p2n14 3,714 s, single run²) | **297×** | 347× | 73.9× | 141× |

- **Correctness:** every COAL config matched truth or the reference on every case run: Tier A 412/412, Tier C
  2000/2000 simulated plus 3273/3273 NCBI cases matching p2n14 (excluding invalid inputs), and **Tier D 2595/2595
  real Compara trees matching truth**. On those trees p2n14 failed 22 (19 selenocysteine, 2 open, 1 regex
  capture-group limit) and p2n16 failed 2 (the same 2 open trees, not investigated). COAL's only output
  difference from v14 is the documented selenocysteine fix.
- **coal-cpp is the fast path everywhere.** coal-py is fast on large inputs, but **slower than p2n14 on small
  files** (Python start-up and imports: 92 ms baseline against 24 ms for p2n14).
- **Thread scaling is weak.** coal-cpp plateaus at about 1.2× from `-j2` (efficiency 6% at -j20). coal-py gains 1.9× at default
  `-j`, but its committed memory goes from 0.69 GB to about 10 GB at n100k (known COAL issue, below).
- At 1M: p2n14 **measured, 3,714 s and 18.6 GB, TRUTH** (one overnight run, no warm-up). p2n16 **DNF (>900 s)**,
  so the speed-ups against it are lower bounds.

² p2n14 at 1M is one run with no warm-up, made overnight (user decision; NOTES.md). The COAL figures are medians of
3 after a warm-up, except coal-py default (single run). The earlier extrapolation of about 2,400 s
(about 190×/220×/48×/90×) **underestimated p2n14's time by about 1.5×** (measured/extrapolated = 1.55). See the 1M
section.

Plots: `tierB_nsweep_time.png`, `tierB_nsweep_mem.png`, `tierB_lsweep_time.png`, `tierE_coalcpp_scaling.png`,
`tierC_ncbi100_per_file.png`. Merged correctness for every tier: `correctness.csv`.

## Tier B: scaling (warm medians, seconds)

| case | p2n14 | p2n16 | coal-py -j1 | coal-py default | coal-cpp -j1 | coal-cpp default |
|---|---|---|---|---|---|---|
| n1000_L400 | 2.31 | 0.434 | 0.185 | 0.180 | 0.031 | 0.027 |
| n10000_L400 | 25.2 | 4.43 | 0.736 | 0.559 | 0.130 | 0.112 |
| n100000_L400 | 239 | 98.9 | 5.04 | 2.64 | 1.07 | 0.885 |
| n1000_L100 | 0.634 | 0.178 | 0.134 | 0.134 | 0.026 | 0.025 |
| n1000_L1000 | 6.46 | 0.915 | 0.272 | 0.265 | 0.048 | 0.039 |
| n1000_L10000 | 280 | 10.4 | 0.958 | 0.861 | 0.251 | 0.216 |
| n10_L35000 | ERROR¹ | 0.628 | 0.157 | 0.158 | 0.029 | 0.026 |
| n100_L35000 | ERROR¹ | 5.65 | 0.604 | 0.646 | 0.125 | 0.086 |

5 timed runs after 1 graded warm-up. The default `-j` is 20 for coal-cpp and automatic for coal-py. -j12 is an
extra, not a headline; it is in `runs.csv` and within 6% of the default everywhere. Drift check (coal-cpp -j1
n10000_L400): 0.1279 s at the start, 0.1295 s at the end. Cold-cache runs (RAMMap standby flush, verified): `runs.csv`
with `cold=1`.

¹ **p2n14 cannot process proteins this long. This is a genuine v14 limit, not a harness or platform issue.** v14
builds one regex per sequence with 1–7 capture groups per residue, and Perl allows at most 65,535 groups per regex.
These proteins need 96,000–102,506. The threshold depends on composition: about 22,000–24,000 aa for typical
proteins. The failure was reproduced outside the harness, identically, under Strawberry Perl 5.42.3, Cygwin Perl
5.42.3 and Linux Perl 5.40.1. Details: NOTES.md.

### 1M sequences (n1000000_L400, 1.2 Gnt)

| tool | result | peak committed memory |
|---|---|---|
| coal-cpp -j1 | 12.52 s (median of 3) | 172 MB |
| coal-cpp default | 10.71 s (median of 3) | 172 MB |
| coal-py -j1 | 50.27 s (median of 3) | 1,176 MB |
| coal-py default | 26.43 s (**single run**, under the 24 GB cap: OK, no OOM) | 8,032 MB |
| p2n16 | **TIMEOUT** (single run, no warm-up, 900 s limit set for the deadline) | 3,743 MB at the kill |
| p2n14 | **3,714.1 s** (62 min; **single run, no warm-up**, overnight, 3 h timeout), **TRUTH** | **18,563 MB** (under the 24 GB cap: OK) |

p2n14 at 1M grows faster than linearly: 15.5× the n100k time for 10× the sequences, an exponent of 1.19 for that
decade against 0.98 for the one before. Memory grows 10.0× per decade (1,852 MB → 18,563 MB). The cause of the
superlinear time was not investigated.

**Extrapolations, kept for comparison. Not measurements.**
- p2n14 at 1M, made before the overnight run: **about 40 min (~2,400 s), about 18–19 GB.** It assumed linear time
  (n-sweep exponents 1.04 and 0.98) and about 9.6× memory per decade. **Time was underestimated by about 1.5×**
  (3,714 s measured). The memory estimate was right (18.6 GB measured).
- **p2n16 at 1M:** measured TIMEOUT at 15 min. Its n-sweep is superlinear (99 s at n100k, exponent 1.35 in the last
  decade), which extrapolates to about 35–40 min (99 s × 10^1.35), consistent with the timeout. Its memory trend (6.7 GB at n100k, about 9.8× per decade)
  suggested an OOM at about 65 GB. It hit the time limit first, at 3.7 GB, so whether it would OOM under the cap is
  unknown.
- coal-py default at 1M: the user expected an OOM and my code-based estimate was about 13–14 GB. **Measured: 8.0 GB,
  OK.** See the memory issue below.

### Tier E (reduced): coal-cpp threads (`tierE_coalcpp_scaling.png`)

Speedup of coal-cpp over `-j1`, 1 warm-up + 3 runs:

| -j | 1 | 2 | 4 | 8 | 12 | 20 |
|---|---|---|---|---|---|---|
| n100000_L400 (s) | 1.195 | 0.989 | 0.976 | 0.957 | 1.004 | 1.018 |
| speed-up | 1.00 | 1.21 | 1.22 | 1.25 | 1.19 | 1.17 |
| n1000000_L400 (s) | 12.42 | 10.17 | 10.74 | 10.85 | 10.88 | 10.83 |
| speed-up | 1.00 | 1.22 | 1.16 | 1.14 | 1.14 | 1.15 |

coal-cpp stops speeding up at `-j2`: parallel efficiency is 61% at -j2 and 6% at -j20. The P-core (6), P-thread (12)
and physical-core (14) marks on the plot make no difference, because the plateau starts long before them. So most
of coal-cpp's time is in a part of the run that does not use extra threads. It is not a hardware limit. Which part
it is was not profiled. coal-py on
n100000_L400: -j1 4.91 s / 689 MB, default 2.56 s / 10.3 GB (1.92×).

### Tier C (reduced): per-file loop on real data (`tierC_ncbi100_per_file.png`)

One fixed-seed batch of 100 NCBI families, drawn from the 3273 cases every config passed (INPUT_INVALID and the
p2n14-failing SELD/FDHF excluded). One process per file, 1 warm-up + 3 timed. Median time for the batch:

| p2n14 | p2n16 | coal-py -j1 | coal-py default | coal-cpp -j1 | coal-cpp default |
|---|---|---|---|---|---|
| 12.52 s (125 ms/file) | 9.05 s (90) | 18.43 s (184) | 18.44 s (184) | 1.84 s (18) | 1.84 s (18) |

Batch time is linear in the number of files at this scale, so the other batch sizes were dropped. Also measured
before the stop: simulated families, all 2,000 files, coal-cpp default 41.45 s (21 ms/file), tagged INCOMPLETE
because the other configs did not run. Native batch mode (mode 2) is N/A for COAL. The in-process API loop was not
run.

### Tier D: Ensembl Compara 116 trees (graded against truth)

**Correctness** (`correctness_D.csv`, `tierD_table.md`). 10 workers; the 10 largest trees run serially. Stopped at the
21:00 cut, so **2,595 of 2,634 trees** were graded. **bin4: 250/289; the remaining 39 trees are not run (time)**,
including 1 of the 10 largest.

| bin | trees | p2n14 | p2n16 | coal-py -j1 / default | coal-cpp -j1 / default |
|---|---|---|---|---|---|
| bin0 2–10 | 482/482 | 481 T, 1 E | 482 T | 482 T | 482 T |
| bin1 11–50 | 449/449 | 449 T | 449 T | 449 T | 449 T |
| bin2 51–200 | 656/656 | 649 T, 7 E | 656 T | 656 T | 656 T |
| bin3 201–1000 | 758/758 | 747 T, 11 E | 757 T, 1 E | 758 T | 758 T |
| bin4 1001+ | 250/289 | 247 T, 3 E | 249 T, 1 E | 250 T | 250 T |
| **all** | 2595/2634 | 2573 T, 22 E | 2593 T, 2 E | **2595 T** | **2595 T** |

**COAL matches truth on every graded tree** (0 failures in 10,380 COAL runs).

**Why PAL2NAL failed** (p2n14 rerun directly to capture stderr):
- **p2n14, 22 failures: 19 selenocysteine + 2 open + 1 capture-group limit.**
  - **19 trees, selenocysteine.** Every one has `U` in its proteins (2–225 per tree; 19 of the 38 selenocysteine
    trees graded). This is the documented v14 defect (`fix_selenocysteine`): "inconsistency", **exit 0**, no output.
    The stop_in_protein / x_residue / ambiguous_bases flags on these trees are incidental; the per-flag counts in
    `correctness_D.csv` overlap.
  - **2 trees, t006204 and t000088: open, not investigated.** "inconsistency", exit 0. No `U` and no internal stop;
    flags ambiguous_bases + x_residue.
  - **1 tree, t054129: Perl's 65,535 capture-group limit** (exit 255), the same v14 limit as Tier B L35000, here
    on a real long protein.
- **p2n16, 2 failures: open, not investigated.** The same two trees, t006204 and t000088 ("inconsistency", exit 1).
  Both PAL2NAL versions fail them and COAL matches truth on both, so the cause is probably shared PAL2NAL logic.
  Not confirmed.

**Timing** (`tierD_timing.md`): 20 fixed-seed trees per bin plus the 10 largest, all 6 configs, 1 warm-up + 3. It ran
21:02–21:30. The table gives the total of per-tree medians, with throughput in Mnt/s of CDS:

| bin (20 trees each) | p2n14 | p2n16 | coal-py -j1 | coal-py default | coal-cpp -j1 | coal-cpp default |
|---|---|---|---|---|---|---|
| bin0 2–10 | 0.40 s (0.12) | 1.32 s (0.04) | 2.03 s (0.02) | 2.03 s (0.02) | **0.19 s (0.25)** | 0.19 s (0.25) |
| bin1 11–50 | 1.39 s (0.38) | 1.52 s (0.35) | 2.47 s (0.21) | 2.44 s (0.22) | **0.34 s (1.53)** | 0.35 s (1.50) |
| bin2 51–200 | 10.14 s (0.42) | 3.06 s (1.40) | 3.76 s (1.14) | 3.79 s (1.13) | **0.44 s (9.62)** | 0.43 s (9.88) |
| bin3 201–1000 | 46.12 s (0.29) | 9.86 s (1.34) | 5.90 s (2.24) | 5.91 s (2.23) | **0.69 s (19.1)** | 0.60 s (22.0) |
| bin4 1001+ | 184.44 s (0.18) | 59.84 s (0.55) | 25.85 s (1.27) | 25.81 s (1.27) | **2.74 s (12.0)** | 2.03 s (16.1) |
| 10 largest | **not run (time)** | not run (time) | not run (time) | not run (time) | not run (time) | not run (time) |

- coal-cpp -j1 is 2–67× faster than p2n14 depending on the bin (bin0 2.1×, bin3 67×). coal-py is **slower** than
  p2n14 on small trees (bin0–bin1, start-up cost) and faster from bin2 on (2.7× in bin2, 7.1× in bin4).
- **10 largest: not run (time).** Only the largest tree (t000052, 7.58 Mnt) was partly timed before the 21:30 stop:
  coal-cpp -j1 0.185 s and default 0.145 s, coal-py -j1 1.65 s. The other three configs were not timed. These
  partial rows are in `runs.csv` but not in the table.

### Scaling observations

- **COAL's thread scaling is modest.** At n100000_L400, the default `-j` over `-j1` is **1.2× for coal-cpp**
  (1.07 → 0.885 s) and **1.9× for coal-py** (5.04 → 2.64 s). On the other Tier B cases coal-cpp gains 1.16–1.45×.
  -j12 (one per P-core thread) is no better than -j20. See Tier E: coal-cpp levels off at `-j2`.
- **coal-py's default `-j` costs memory (known COAL issue; GitHub issue draft in NOTES.md).** At n100000_L400
  its peak committed memory (summed over the job) is 9.0–10.8 GB at the default `-j` against 0.69 GB at `-j1`:
  about 15× the memory for 1.9× the speed. At 1M it is 8.0 GB. From reading the code: each process imports numpy
  once it has checked 2 M codons (`core.py:353-371`), and that adds hundreds of MB per process. The exact cost per
  process, and whether numpy's OpenBLAS thread buffers cause it, are unverified (the 1M result does not fit a fixed
  cost). Not fixed in this run.
- **coal-py's default `-j` can be slightly slower than `-j1`.** On n100_L35000 it is 0.94× (automatic `-j`
  resolves to 2 workers there).
- **p2n16 scales superlinearly in N.** Time per 1,000 sequences is 434 ms at n1k, 443 ms at n10k and 989 ms at
  n100k (per-decade exponent 1.01, then 1.35). Peak memory reaches 6.7 GB at n100k. The cause was not
  investigated. p2n14 is linear in N (exponents 1.04, 0.98) but superlinear in length (exponent 1.64 from L1000 to
  L10000). COAL is sublinear over these ranges because start-up dominates the small cases.

## Correctness

### Tier A
412/412 TRUTH for every config. p2n14 prints spurious warnings on table 10 but its output is correct
(NOTES.md).

### Tier C grading (before timing)

| set | p2n14 | p2n16 | coal-py (-j1, default) | coal-cpp (-j1, default) |
|---|---|---|---|---|
| families (2000, simulated, with truth) | 2000 TRUTH | 2000 TRUTH | 2000 TRUTH | 2000 TRUTH |
| ncbi (3437, real, reference = p2n14) | 3273 REF, 2 ERROR, 162 INPUT_INVALID | 3273 REF, 2 KNOWN_DIV, 162 INPUT_INVALID | same as p2n16 | same as p2n16 |

**INPUT_INVALID: 162 of 3437 NCBI cases (4.7%) are bad inputs. This is a test-set construction issue to fix
upstream.** Every tool rejects them: COAL and p2n16 exit 1, while p2n14 exits 0 with empty output (graded ERROR,
as for any tool). They are excluded from timing for every tool. Per-case reasons are in `input_invalid.csv`
(produced by `benchmarks/triage_inputs.py`). The data was not edited.

- **105 cases: duplicate sequence IDs.** IDs are `organism|protein accession`. Bacterial `WP_` accessions are
  shared across strains, and the organism name is generic (99 cases "Wolbachia endosymbiont", 5 "Spiroplasma
  endosymbiont", 1 "Rickettsiella endosymbiont"). So different genomes carrying the same protein get the same ID
  (e.g. RPMB: 233 sequences from 233 distinct genomes, 2 colliding IDs). COAL rejects the file ("duplicate sequence
  ID"). p2n16 merges the duplicates into double-length rows ("rows are not all the same length"). Upstream fix: add
  the genome/assembly accession to the ID, or deduplicate identical proteins. Scenarios: table11_large 92,
  table11_medium 8, table4_medium 5.
- **57 cases: the protein is not the translation of its CDS** for at least one sequence (residues the CDS does not
  encode, or the CDS is shorter than 3× the protein). All 57 carry an NCBI `exception` and/or
  `cds_len_not_3x_protein` flag, so the builder could filter them out or move them to a negative-test set.
  Scenarios: table1_small 44, table1_medium 12, table11_large 1.

### Input validation (how each tool reports a bad input)

On the 105 duplicate-ID inputs:
- **COAL (coal-py and coal-cpp) reports it clearly.** It exits 1 with
  `coal: error: <pep file>: duplicate sequence ID '<id>'`, naming the file and the ID.
- **p2n16's error is misleading.** It exits 1 with "the alignment rows are not all the same length; most have N
  columns", then lists the duplicated ID with a doubled length. It merges duplicate IDs silently before checking, so
  the error points at row length, not at the real cause.
- **p2n14 fails silently, as far as the exit code shows.** It exits **0** with empty output. A pipeline that checks
  exit codes would carry on with an empty alignment.

On the 57 protein ≠ CDS inputs, COAL names the sequence and the column where the DNA stops matching. p2n16 prints
its "inconsistency" block and exits 1. p2n14 exits 0 with empty output (its stderr on these was not kept; on SELD/FDHF
it is the same "inconsistency" block).

**KNOWN_DIV: SELD and FDHF (selenocysteine, table 11): v14 fails; the others agree.** v14 accepts `U` but
gives it no codon, prints "inconsistency between pep and nuc", exits 0 and writes nothing (documented: p2n16
`tests/divergences.tsv` `fix_selenocysteine`; COAL README "U (selenocysteine)"). p2n16, coal-py and coal-cpp (both
`-j` settings) produce byte-identical output on both cases, with TGA at all 29 selenocysteine positions. Because
p2n14 has no output there, these two cases are not timed for any tool.

## Method notes and caveats

- **The 24 GB memory cap is artificial.** Every run is in a Windows Job Object limited to 24 GB of committed memory
  (the machine has 32 GB). A run that hits it would be recorded as `OOM (cap=24GB)`. That would mean "needs more
  than 24 GB here", not a crash on a bigger machine. No run hit the cap. Memory is peak committed memory summed over
  the job (all child processes), which can be larger than resident memory.
- Timing: warm medians of timed runs after a graded warm-up. Every timed output is checked byte-for-byte against the
  graded one, and a (tool, case) that failed grading is not timed. Runs are strictly sequential with nothing else
  timed running. Clock speed is logged (130–135% of nominal throughout); the coal-cpp drift check stayed within
  0.126–0.138 s across all steps.
- Cold-cache runs (Tier B, RAMMap standby flush, verified before and after) are in `runs.csv` with `cold=1`.
- `harness_commit` is recorded per row, separately from `tool_commit`. See NOTES.md for history and the false
  `-dirty` label on some rows.
- Defender exclusion on `C:\coalbench`, High-performance power plan. Background load is listed in NOTES.md.

## Scope reductions due to time

A hard deadline came in at 19:5x (results and PR by 21:45, later moved to 22:00). These cuts were made on the user's instruction:

1. **Runs per measurement: 1 warm-up + 3 timed** for everything after 19:57 (1M, Tier E, reduced Tier C, Tier D),
   instead of the planned 5 (Tier B before the deadline used 5).
2. **Tier C timing stopped at 19:57.** The full design (batch sizes 1/10/100/1000/all on both sets, 5 configs) was
   replaced by one fixed-seed batch of 100 NCBI families, all 6 configs. The families/2000 batch measured before
   the stop (coal-cpp default only) is kept and tagged INCOMPLETE. Batch sizes 1, 10, 1000 and all: **not run**. The
   in-process API loop: **not run**.
3. **1M: COAL only.** coal-cpp -j1/default and coal-py -j1, timed. coal-py default: one run. p2n16: one run, no
   warm-up, 900 s timeout (TIMEOUT). p2n14 at 1M: not run before the deadline. It was **run afterwards, overnight**
   (one run, no warm-up; see the 1M section). Cold-cache 1M: **not run**.
   -j12 at 1M: only in Tier E (coal-cpp).
4. **Tier E reduced** to coal-cpp at -j 1, 2, 4, 8, 12, 20 on n100000_L400 and n1000000_L400, plus coal-py -j1 and
   default on n100000_L400. -j 6, 10, 14, 16, the 10 largest Compara trees and the 1,000-family batch: **not run**.
5. **Tier D correctness** was cut at 21:00: 2,595/2,634 trees graded. **bin4: 39 trees not run (time).**
6. **Tier D timing** was reduced to 20 trees per bin plus the 10 largest (instead of up to 50 per bin plus the 10
   largest), 1 warm-up + 3. It ran 21:02–21:30 and stopped at 21:30 wherever it had got to. Bins not reached are
   **not run (time)** (see Tier D).
7. The deadline moved from 21:45 to 22:00 at about 20:45, and the stop point moved from 21:15 to 21:30.
8. Planned but **not run**: the 1M case in the cold-cache block, Tier B -j12 at 1M (except coal-cpp in Tier E), and a
   coal-py thread curve.
