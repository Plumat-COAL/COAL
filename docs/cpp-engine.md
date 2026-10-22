# C++ engine (coal-cpp)

`coal-cpp` is a C++17 implementation of the same engine (standard library only, MIT). It takes the same options as
`coal` and writes the same output, warnings, report and exit codes, so it can replace `coal` in any command line.
`coal` itself hands inputs of 1 MB or more to `coal-cpp` when it finds the same version (set `COAL_CPP` to its path,
or put it next to Python or on `PATH`; `COAL_ENGINE=python|cpp` chooses), and the Windows app ships it. To build it:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
build/coal-cpp pep.aln nuc.fasta -o codon_aln.fasta       # Windows (Visual Studio): build\Release\coal-cpp.exe
```

- **Same behaviour as the Python engine:** all 285 stored PAL2NAL v14 reference checks pass with `COAL_IMPL=cpp`;
  3,000 random inputs (all formats and options, report JSON included) give identical results in both engines for
  two seeds; and the random comparison with the real PAL2NAL v14 (`tests/fuzz_vs_pal2nal.py`, 1,500 cases per seed,
  seeds 20261017, 777 and 4242) finds 0 unexplained differences (these larger runs were done locally). CI builds it
  with warnings as errors on Linux, Windows and macOS and runs the reference checks and 1,000 random C++-vs-Python
  inputs (one seed) on every push.
- **Speed:** 3-6x faster than the Python engine, with a third to a half of its peak memory (measured on a 2-core
  Linux VM, commit `4084d1f`; `benchmarks/results/cpp_vs_python_linux_4084d1f.md`):

| Sequences × codons | Python | C++ | C++ speed-up | Peak memory Python / C++ |
|---|---:|---:|---:|---:|
| 10,000 × 500 | 0.651 s | 0.109 s | 6.0× | 57 / 20 MB |
| 100 × 10,000 | 0.175 s | 0.028 s | 6.3× | 27 / 14 MB |
| 100,000 × 500 | 2.945 s | 0.964 s | 3.1× | 125 / 37 MB |

  (Default `-j` for both engines; outputs byte-identical. `--impl both` in `benchmarks/run_bench.py` reproduces it.)
- **Tests:** `COAL_IMPL=cpp COAL_CPP_BIN=build/coal-cpp python -m unittest tests.test_differential
  tests.test_cpp_vs_python tests.test_cpp_cli` (the C++ tests skip when the binary isn't built).
