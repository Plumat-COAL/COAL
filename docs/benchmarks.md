# Benchmark results

## COAL 47f57f1 vs PAL2NAL v14 (Windows, up to 1 million sequences)

Speed-up of COAL over PAL2NAL v14, from warm medians (5 timed runs after a warm-up; 3 for the 1M case):

| Workload | PAL2NAL v14 | coal-cpp, 1 thread | Speed-up | coal-py, 1 process | Speed-up |
|---|---:|---:|---:|---:|---:|
| 1,000,000 sequences × 400 codons | 3,714 s (62 min), 18.6 GB | 12.5 s, 172 MB | **297×** | 50.3 s, 1.2 GB | 74× |
| 100,000 sequences × 400 codons | 239 s | 1.07 s | **223×** | 5.04 s | 47× |
| 1,000 sequences × 10,000 codons | 280 s | 0.251 s | **1,115×** | 0.958 s | 292× |
| 1,000 sequences × 400 codons | 2.31 s | 0.031 s | 74× | 0.185 s | 12.5× |
| 100 sequences × 35,000 codons | fails¹ | 0.125 s | | 0.604 s | |

- Every COAL output matched the simulated truth or PAL2NAL v14, here and on 2,595 real Ensembl Compara gene trees
  (where PAL2NAL v14 failed on 22).
- The PAL2NAL v14 run at 1M is a single run with no warm-up (62 minutes). With all threads (default `-j`), coal-cpp
  is 347× faster there; extra threads help it little (about 1.2× from 2 threads on).
- On many small files, one process per file, the Python engine was slower than PAL2NAL v14 (184 against 125 ms per
  file, Python start-up); coal-cpp took 18 ms per file.
- ¹ PAL2NAL v14 builds one regex per sequence and Perl allows at most 65,535 capture groups in one, so it fails on
  proteins longer than about 22,000–24,000 residues (reproduced on three Perl builds).

Machine: Intel i5-13600KF (14 cores / 20 threads), 32 GB, Windows 11. Memory is peak committed memory. Full results,
including PAL2NAL v16 and the method: `benchmarks/results/windows_47f57f1/SUMMARY.md` on the `bench/windows-47f57f1`
branch.

The sections below are from earlier versions of the Python engine (the commit is named in each section).

## COAL vs PAL2NAL v14 (v0.3)

COAL vs PAL2NAL v14 on synthetic datasets (see `benchmarks/README.md`). Each run is one
command-line process, timed from start to finish (including Perl/Python start-up).
Times are the median of 3 runs; COAL ran as a single process (`-j 1`), with the Python engine at an earlier
commit (`58a1219`, see below). **On these benchmark datasets every COAL output was byte-identical to PAL2NAL v14's.**

| Sequences × codons | Total codons | PAL2NAL | COAL | Speedup | Peak memory PAL2NAL / COAL |
|---|---:|---:|---:|---:|---:|
| 10 × 500 | 5 k | 0.051 s | 0.061 s | 0.8× | 13 / 15 MB |
| 100 × 100 | 10 k | 0.077 s | 0.062 s | 1.3× | 11 / 15 MB |
| 100 × 500 | 50 k | 0.265 s | 0.072 s | 3.7× | 15 / 15 MB |
| 100 × 2,000 | 200 k | 1.36 s | 0.093 s | 14.6× | 23 / 17 MB |
| 1,000 × 500 | 500 k | 2.39 s | 0.139 s | 17.3× | 35 / 19 MB |
| 100 × 10,000 | 1 M | 24.5 s | 0.207 s | **118.6×** | 55 / 23 MB |
| 10,000 × 500 | 5 M | 24.7 s | 0.864 s | **28.6×** | 219 / 55 MB |

What the numbers show:

- **More sequences:** COAL is 17–29× faster from 1,000 sequences upward.
- **Longer genes:** PAL2NAL's time grows much faster than the gene length: 5× more codons
  (2,000 → 10,000) took 18× longer. COAL's time grows in proportion to the input, so the gap
  widens with gene length.
- **Tiny inputs:** both tools take about 0.05 s, which is mostly start-up time.

Machine: Windows 11, Intel 64-bit CPU (Family 6 Model 183), Python 3.13.9, Strawberry Perl 5.42.3, COAL commit
`58a1219`. Raw data: `benchmarks/results/bench_windows_20261016_2356.{csv,txt}`.

**Limitations:** synthetic data only (random ancestral gene, about 10% of codons changed, about 5% deleted; no UTRs or
frameshifts in the benchmark data, although COAL handles both); one machine; 3 repeats. Speedups depend on the
machine and dataset.

## Scaling (v0.3)

Synthetic data, 500 codons per gene, on the same Windows 11 machine as above (20 logical CPUs). One run per
configuration. All outputs were byte-identical.

| Sequences | Input size | 1 process | 4 processes | All cores (20) | Peak memory, 1 process |
|---:|---:|---:|---:|---:|---:|
| 100,000 | 194 MB | 7.1 s | 2.4 s | 2.5 s | 55 MB |
| 500,000 | 974 MB | 25.2 s | 9.8 s | 10.0 s | 131 MB |

- **Memory** stays small: v0.2 needed about 2.5 GB for the 500,000-sequence input (measured on Linux).
- **Throughput** is about 100 MB of input per second with 4 processes. PAL2NAL was not run at these sizes; from its
  10,000-sequence time (24.7 s), and assuming it grows in proportion, it would need roughly 4 and 21 minutes.
- **More than 4 processes don't help yet.** The main process reads all the DNA, passes it to the workers and writes
  all the output, and that is now the bottleneck.

Raw data: `benchmarks/results/bench_windows_20261017_0059.csv` (500,000) and the 100,000 row of
`bench_windows_20261017_0055.csv`. That file's 500,000 row is not used: its 1-process column was accidentally
timed with automatic `-j`. Memory is measured for the main process only.

On the 2-core Linux development machine (v0.2 → v0.3, 1 process): 100,000 sequences 6.3 s → 3.6 s,
515 → 87 MB; 500,000 sequences 32.0 s → 17.0 s, 2,523 → 143 MB.

Reproduce with `python benchmarks/run_bench.py --no-pal2nal --sizes 100000x500 500000x500 --coal-jobs 4 0`.

## Reproduce

```
python scripts/fetch_pal2nal.py            # once; PAL2NAL is not redistributed
pip install -e . psutil                    # psutil only for the memory columns
python benchmarks/run_bench.py --repeats 3 --perl /path/to/perl
```

The script generates the same datasets every time (fixed random seed), times both tools,
checks their outputs match, and writes `benchmarks/results/<name>.csv` plus a `.txt` file
with the machine and software versions.
