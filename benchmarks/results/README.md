# Benchmark results

<!-- Superseded runs are in the git history. Every file left here is cited by README.md or pitch/PITCH.md. -->
Each run is a `.csv` (one row per dataset) and a `.txt` (machine, versions). Only compare numbers from the same machine.

| Files | What | Cited in |
|---|---|---|
| `bench_windows_20261016_2356.*` | COAL vs PAL2NAL v14, Windows 11, 20 logical CPUs, COAL `58a1219`: up to 119x faster, outputs identical | docs/benchmarks.md |
| `bench_windows_20261017_0055.*`, `bench_windows_20261017_0059.*` | COAL scaling to 100,000 and 500,000 sequences, 1 / 4 / 20 processes (same machine) | README "Scaling (v0.3)" |
| `windows_47f57f1/` | COAL vs PAL2NAL v14 and v16, Windows 11, i5-13600KF, COAL `47f57f1`: the README's 297× (1 million sequences) and 1,115× (long genes) figures; `SUMMARY.md` first. `correctness.csv.gz` is the merged per-case correctness table; the per-tier copies and the 20 per-gene `diffs/` (17 MB) stay on the `bench-windows-47f57f1` tag | README, docs/benchmarks.md |
| `cpp_vs_python_linux_4084d1f.md` + `cpp_vs_python_linux_4084d1f_*` | Python vs C++ engine, 2-core Linux VM, commit `4084d1f`: C++ 3-6x faster, a third of the memory | docs/cpp-engine.md |

Reproduce with `benchmarks/run_bench.py` (see `benchmarks/README.md`).
