# Protein-alignment layer benchmark

Machine: Linux-6.18.44-fc-v378-x86_64-with-glibc2.39, 2 hardware threads, Python 3.13.18. DIAMOND: diamond version 2.2.8.
Data: 1000 + 1000 CDS -> 999 DIAMOND reciprocal best hits; 313.8 M DP cells; median protein length 384.
DIAMOND hit coverage of the shorter protein: median 1.000, 35 pairs < 0.9, 6 pairs < 0.5.

## Time (protein alignment of all pairs; median of 5 runs)

| Method | Threads | Wall (s) | CPU (s) | Pairs/s (wall) | GCUPS (CPU) |
|---|---:|---:|---:|---:|---:|
| coal-cpp-align NW, needle defaults | 1 | 0.643 | 0.426 | 1,553 | 0.74 |
| coal-cpp-align NW, needle defaults | 2 | 0.644 | 0.410 | 1,550 | 0.76 |
| coal-cpp-align NW, pwalign preset | 1 | 0.678 | 0.473 | 1,473 | 0.66 |
| coal-cpp-align NW, build 'baseline-x86-64' | 1 | 1.000 | 0.657 | 999 | 0.48 |
| coal-cpp-align NW, build 'no-simd' | 1 | 2.170 | 1.439 | 460 | 0.22 |
| coal-cpp-align NW, build 'clang-avx2' | 1 | 0.477 | 0.318 | 2,094 | 0.99 |
| Biopython PairwiseAligner (same scoring) | 1 | 5.849 | 3.969 | 171 | 0.08 |
| parasail sg_trace_striped_16 (open 10, extend 1) | 1 | 1.323 | 0.884 | 755 | 0.36 |
| parasail sg_trace_scan_16 (open 10, extend 1) | 1 | 0.322 | 0.211 | 3,100 | 1.48 |
| coal-cpp-align NW, DIAMOND's gaps (open 11, extend 1, free ends) | 1 | 0.598 | 0.391 | 1,671 | 0.80 |
| clustalo 1.2.4, one process per pair | 2 | 37.957 | 30.926 | 26 | 0.01 |
| mafft v7.526 (2024/Apr/26), one process per pair | 2 | 179.815 | 205.577 | 6 | 0.00 |
| DIAMOND hit expanded, flanks=gap | 1 | 0.036 | 0.023 | 27,697 | 13.81 |
| DIAMOND hit expanded, flanks=align | 1 | 0.050 | 0.032 | 20,155 | 9.71 |

DIAMOND blastp both ways (2 threads), RBH input: 2.04 s with `qseqid sseqid bitscore`, 2.22 s with the gapped columns added.

## dN/dS (NG86, Biopython) against the needle-default NW

Codon alignments by coal-cpp; codons with a gap or stop in either sequence dropped. r = Pearson over pairs with finite values in both; |dw| = |dN/dS - reference dN/dS|.

| Method | codons compared (sum) | r dN | r dS | r dN/dS | median abs dw | pairs abs dw > 0.05 | NA dS |
|---|---:|---:|---:|---:|---:|---:|---:|
| NW needle (reference) | 445,281 | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0 | 1 |
| NW pwalign preset | 445,838 | 0.9545 | 0.9433 | 0.9920 | 0.0000 | 33 | 1 |
| NW open 11 / extend 1 | 445,524 | 0.9905 | 0.9925 | 0.9978 | 0.0000 | 13 | 1 |
| clustalo (pair as a cluster) | 445,707 | 0.9767 | 0.9873 | 0.9883 | 0.0000 | 77 | 1 |
| mafft (pair as a cluster) | 445,573 | 0.9898 | 0.9932 | 0.9962 | 0.0000 | 23 | 1 |
| DIAMOND gapped, flanks=gap | 440,800 | 0.8837 | 0.8913 | 0.9359 | 0.0000 | 99 | 1 |
| DIAMOND gapped, flanks=align | 445,099 | 0.9321 | 0.9341 | 0.9715 | 0.0000 | 81 | 1 |

By DIAMOND hit coverage of the shorter protein (r dN/dS against the reference; pairs):

| Method | coverage >= 0.9 | coverage < 0.9 |
|---|---:|---:|
| NW needle (reference) | 1.0000 (964) | 1.0000 (34) |
| NW pwalign preset | 0.9948 (964) | 0.9613 (34) |
| NW open 11 / extend 1 | 0.9985 (964) | 0.9917 (34) |
| clustalo (pair as a cluster) | 0.9902 (964) | 0.9683 (34) |
| mafft (pair as a cluster) | 0.9972 (964) | 0.9870 (34) |
| DIAMOND gapped, flanks=gap | 0.9711 (964) | 0.6000 (34) |
| DIAMOND gapped, flanks=align | 0.9728 (964) | 0.9503 (34) |

## Supplementary: DIAMOND with `--masking 0` (same data, same NG86 comparison)

DIAMOND masks low-complexity regions by default, and inside a masked stretch its alignment is arbitrary: of the 964
pairs whose hit covers >= 90% of the shorter protein, 409 have a block that differs from the optimal global alignment
under DIAMOND's own gap costs (open 11, extend 1) over the same region, 219 of them with masked residues. Rerunning
blastp with `--masking 0` (999 of 999 RBH pairs still hit):

| Method | r dN/dS | pairs abs dw > 0.05 |
|---|---:|---:|
| DIAMOND gapped, `--masking 0`, flanks=gap | 0.9825 | 53 |
| DIAMOND gapped, `--masking 0`, flanks=align | 0.9881 | 45 |

Reproduce: `python benchmarks/align_bench.py ... --masking 0`.

Notes: run 2026-11-11 in a 2-core cloud VM shared with other jobs, so wall times are noisy and the 2-thread rows
show no scaling; CPU time (last-but-one column) is the steadier number. Builds: default = GCC 13 -O3 with the AVX2
row loops picked at run time; 'baseline-x86-64' = -DCOAL_ALIGN_SCALAR (SSE2 vectorised loops, scalar scan; what
MSVC and non-AVX2 CPUs get); 'no-simd' = also -fno-tree-vectorize; 'clang-avx2' = Clang 18. parasail 2.6.1 (Python
wheel) uses 16-bit lanes and integer gaps (open 10, extend 1 here), so it is a speed reference, not the same scoring.
