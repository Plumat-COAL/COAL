# COAL engines after the merge: Python vs C++ (measured, run 5 step 4)

Commit **4084d1f** (`overnight/core`: Keanu's merge 0729829 + CPP-25..CPP-33). Produced with
`python benchmarks/run_bench.py --impl both --coal-jobs 0 --sizes <N>x<L> --repeats <R> --data-dir <local disk>`
(one call per size; raw rows in `cpp_vs_python_linux_4084d1f_*.csv`, run notes in the matching `*.txt`).

**Machine:** Linux VM on an Intel Core i5-13600KF host: **2 CPU cores, 4 GB RAM** visible to the VM, Ubuntu 22.04,
Python 3.10.12 (numpy available), g++ 11.4 Release (-O3). Datasets on the VM's local disk. Wall-clock includes
process start-up; medians of 5 runs (3 for 100,000 x 500). Output to stdout (fasta). All four outputs byte-identical
in every size.

Columns: Python -j 1 = `python -m coal -j 1`; Python default = automatic -j (all cores for large inputs);
C++ -j 1 = `coal-cpp -j 1`; C++ default = all hardware threads (2 here).

| dataset (seqs x codons) | DNA MB | Python -j 1 (s) | Python default (s) | C++ -j 1 (s) | C++ default (s) | speed-up -j 1 | speed-up default | peak MB (Py -j1 / Py default / C++ -j1 / C++ default) | identical |
|---|---|---|---|---|---|---|---|---|---|
| 10,000 x 500 | 13.69 | 0.682 | 0.651 | 0.125 | 0.109 | 5.5x | 6.0x | 51.7 / 57.3 / 17.2 / 19.6 | True |
| 100 x 10,000 | 2.72 | 0.176 | 0.175 | 0.033 | 0.028 | 5.3x | 6.3x | 27.7 / 26.6 / 13.8 / 14.2 | True |
| 100,000 x 500 | 136.94 | 4.213 | 2.945 | 1.050 | 0.964 | 4.0x | 3.1x | 105.2 / 125.0 / 35.2 / 36.5 | True |

Notes:
- Only 2 cores: the default -j uses at most 2 workers (Python uses them only for the large input: 4.2 s -> 2.9 s;
  C++ gains 8-15%). On a machine with more cores the default columns will scale further: re-run the same
  command there for headline numbers.
- The first step-4 run (files `cpp_vs_python_linux_47d2168_*`, not committed) found the merged report's
  end-of-run step holding a view of every entry: C++ peak 152 MB on 100,000 x 500. Fixed in CPP-32/CPP-33
  (now 35-37 MB; guarded by tests/test_cpp_vs_python.py TestCppMemory). Times were unchanged by the fix.
- Pre-merge numbers (01:33, `cpp_vs_python_linux.md`, removed in IMP-M-4; see the git history): C++ 100,000 x 500 0.95 s / 26 MB. The merged engine's
  report keeps the DNA IDs as well (as the Python engine does), hence ~10 MB more.
