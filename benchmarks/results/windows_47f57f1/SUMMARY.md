# COAL Windows benchmark: summary (DRAFT, Tiers C–E still to come)

Tools: COAL `47f57f1` (coal-py, coal-cpp), PAL2NAL v14 (Perl, "p2n14"), PAL2NAL v16 grp-bork Python port
`da572f7` ("p2n16"). Machine and toolchain: `env.json`. Raw data: `runs.csv`, `correctness_*.csv`. Run log and
design notes: `status.log`, `NOTES.md`.

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

### Scaling observations

- **COAL's thread scaling is modest.** At n100000_L400, the default `-j` over `-j1` is **1.2× for coal-cpp**
  (1.07 → 0.885 s) and **1.9× for coal-py** (5.04 → 2.64 s). On the other Tier B cases coal-cpp gains 1.16–1.45×.
  -j12 (one per P-core thread) is no better than -j20. Where the time goes (parallel vs serial parts) was not
  profiled; Tier E will give the full -j curve.
- **coal-py's default `-j` costs memory.** At n100000_L400 its peak committed memory (summed over the job) is
  10.8 GB at the default `-j` against 0.69 GB at `-j1`: about 16× the memory for 1.9× the speed.
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

**KNOWN_DIV: SELD and FDHF (selenocysteine, table 11): v14 fails; the others agree.** v14 accepts `U` but
gives it no codon, prints "inconsistency between pep and nuc", exits 0 and writes nothing (documented: p2n16
`tests/divergences.tsv` `fix_selenocysteine`; COAL README "U (selenocysteine)"). p2n16, coal-py and coal-cpp (both
`-j` settings) produce byte-identical output on both cases, with TGA at all 29 selenocysteine positions. Because
p2n14 has no output there, these two cases are not timed for any tool.

<!-- TODO: Tier C timing, Tier D, Tier E, plots, headline speedups, the 24 GB cap note -->
