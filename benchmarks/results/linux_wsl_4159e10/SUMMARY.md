# COAL 4159e10 vs PAL2NAL v14 on Linux (WSL 2): summary

Made during the full verification of 2026-10-24 (branch `verify/full-2026-10-24`, phase 11 and its second pass).
Tools: COAL `4159e10` (the Python engine and `coal-cpp` built with g++ 15.2, `-O2`, the engines are unchanged in the
commits that follow it), PAL2NAL v14 (Perl 5.40). Produced with `benchmarks/run_bench.py` (one call per file; the
command line is in each `.txt`).

**Machine:** the same desktop as `windows_47f57f1` (Intel Core i5-13600KF: 6 P-cores with SMT + 8 E-cores), here as
a WSL 2 guest: Ubuntu 26.04, kernel 6.6.114, **16 logical CPUs and 29 GB** given to the guest, Python 3.13.13,
`coal-cpp` default `-j` = 16 threads.

**Two storage settings.** The files whose name starts with `small`, `100k` and `small_pal2nal` had the data on the
Windows drive mounted into the guest (DrvFs, `/mnt/d`), where file I/O is slow: kept to show the effect. The `ext4_*`
and `1m_*` files had the data on the guest's own disk (ext4), which is the comparable setting; use those.

## Engines, ext4 (`--impl both`, medians of 3 timed runs after a warm-up; 1M: one run)

| Workload | Python, 1 process | Python, all cores | coal-cpp, 1 thread | coal-cpp, all threads | identical |
|---|---:|---:|---:|---:|---|
| 10,000 × 500 | 0.595 s | 0.440 s | 0.108 s | 0.078 s | yes |
| 100 × 10,000 | 0.210 s | 0.231 s | 0.031 s | 0.026 s | yes |
| 100,000 × 500 | 4.18 s | 2.19 s | 0.898 s | 0.749 s | yes |
| 1,000,000 × 400 | 36.7 s | 17.8 s | 8.44 s | 7.53 s | yes |

DrvFs for comparison: 100,000 × 500 took 23.6 s (Python) and 8.43 s (coal-cpp), file I/O bound.

## Against PAL2NAL v14, ext4 (peak memory is summed PSS of the process tree)

| Workload | PAL2NAL v14 | Python engine, 1 process | Speed-up | coal-cpp, 1 thread (from the table above) | Speed-up |
|---|---:|---:|---:|---:|---:|
| 1,000,000 × 400 | 1,035 s (17 min), 14.9 GB | 35.9 s, 345 MB | 28.8× | 8.44 s | 123× |
| 10,000 × 500 | 11.3 s, 192 MB | 0.581 s, 44 MB | 19.5× | 0.108 s | 105× |
| 100 × 10,000 | 12.8 s, 43 MB | 0.215 s, 21 MB | 59.5× | 0.031 s | 413× |

Every PAL2NAL and COAL output pair was identical (`outputs_identical` in the CSVs).

## Against the Windows table (`windows_47f57f1`)

COAL's shape holds: coal-cpp 5–7× faster than the Python engine, all cores helping the Python engine on many
sequences (4.2 s → 2.2 s; 36.7 s → 17.8 s) and not on long genes, memory flat (345 MB at 1M against 1.2 GB there,
after the 2026-10-23 memory work). PAL2NAL v14, however, runs about 3.6× faster under this Linux Perl than it did on
Windows (1,035 s against 3,714 s at 1M; 11.3 s against 24.7 s at 10,000 × 500 in the v0.3 table), so the speed-ups
here are lower than the headline Windows figures: about 120× (coal-cpp) and 29× (Python) at 1M instead of 297× and 74×.
Both are real; which one a user sees depends on how fast Perl is on their machine.

## Files

- `ext4_both.csv` / `.txt`, `1m_both.csv` / `.txt`: both engines, 1 and all cores.
- `ext4_pal2nal.csv` / `.txt`, `1m_pal2nal.csv` / `.txt`: PAL2NAL v14 against the Python engine (1M: 1 process only).
- `small.csv`, `100k.csv`, `small_pal2nal.csv` (and `.txt`): the same on DrvFs.
- Memory at 500,000 × 150 codons in every output format (both engines, `tests/tools/stress_run.py`) is in the
  verification report on the `verify/full-2026-10-24` branch (`verify/logs/p11-memory.jsonl`): Python 190–197 MB,
  coal-cpp 112–149 MB, identical outputs.
