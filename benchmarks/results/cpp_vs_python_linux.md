# COAL engines: Python vs C++ (measured)

Produced by `python benchmarks/run_bench.py --impl both --sizes <N>x<L> --data-dir <local disk>`
(one call per size; raw rows in the matching `cpp_vs_python_linux_*.csv`, machine notes in `*.txt`).

Machine: Linux VM (Ubuntu 22.04, Python 3.10, g++ 11.4 -O3 Release), **2 CPU cores, 4 GB RAM**; datasets on the VM's
local disk (on the Windows-mounted folder every run is ~0.1 s slower because of file I/O).
Wall-clock includes process start-up. Python = `python -m coal -j 1` (uses numpy when available);
C++ 1 thread = `coal-cpp -j 1`; C++ default = all hardware threads (2 here). All outputs byte-identical.

| dataset (seqs x codons) | DNA MB | Python -j 1 (s) | C++ 1 thread (s) | C++ 2 threads (s) | speed-up (C++ default vs Python) | peak MB Python / C++ |
|---|---|---|---|---|---|---|
| 10,000 x 500 | 13.69 | 0.5926 | 0.1142 | 0.0972 | 6.1x | 44.3 / 15.4 |
| 100 x 10,000 | 2.72 | 0.1491 | 0.0272 | 0.0232 | 6.43x | 23.0 / 11.6 |
| 100,000 x 500 | 136.94 | 3.5588 | 0.9522 | 0.8809 | 4.04x | 57.6 / 25.9 |

Medians of 5 runs (3 for 100,000 x 500). C++ peak memory stays under 30 MB on a 137 MB input (streaming: only IDs and file offsets are kept).
