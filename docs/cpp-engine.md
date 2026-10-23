# C++ engine (coal-cpp)

`coal-cpp` is a C++17 implementation of the same engine (standard library only, MIT). It takes the same options as
`coal` and writes the same output, warnings, report and exit codes, so it can replace `coal` in any command line.
`coal` itself hands inputs of 1 MB or more to `coal-cpp` when it finds the same version (`COAL_ENGINE=python|cpp`
forces one engine). The platform wheels of `coal-align` (Linux x86_64 and aarch64, macOS, Windows x64;
[install.md](install.md)) and the Windows app ship it, so most installs have it without a build. To build it
yourself:

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
- **Not in coal-cpp:** `--batch` (many small runs in one process, [usage.md](usage.md#batch-mode)) is a Python-side
  driver; `coal-cpp --batch` says so and exits 2. coal-cpp starts in a few milliseconds, so one process per job is
  cheap, and `coal --batch` never hands a batch to it.
- **Tests:** `COAL_IMPL=cpp COAL_CPP_BIN=build/coal-cpp python -m unittest tests.test_differential
  tests.test_cpp_vs_python tests.test_cpp_cli` (the C++ tests skip when the binary isn't built).

## Where coal looks for coal-cpp

`coal.dispatch.find_cpp()` takes the first of these; `coal --cpp-path` prints the result (and exits 1, saying why,
when there is none or it is another version), so a pipeline can call the binary directly:
`cpp=$(coal --cpp-path) && "$cpp" pep.aln nuc.fasta -o out.fasta`.

1. **`COAL_CPP`**, if set: that path and nothing else (a build of your own, a shared install).
2. **Inside the installed package:** `coal/bin/coal-cpp` (`coal-cpp.exe` on Windows), where the platform wheels
   put it (`packaging/wheels/build_platform_wheel.py`); used only if it is executable. The pure wheel, the sdist and
   a source checkout have no `coal/bin`.
3. **Next to the Python running `coal`:** the interpreter's folder, its `_internal` (the frozen Windows app) and the
   environment's `bin`/`Scripts` (where Bioconda installs it).
4. **`PATH`.**

A found `coal-cpp` is used only when its `--version` is this package's version; otherwise `coal` says so once
(`coal: notice: ... using the Python engine`) and goes on in Python, unless `COAL_ENGINE=cpp` demands it (exit 1).

### How the platform wheels are made

`release.yml` (job `cpp-wheels`) builds `coal-cpp` per platform and `packaging/wheels/build_platform_wheel.py` copies
it into the pure wheel as `coal/bin/coal-cpp` with the executable bit set in the zip, rewrites `RECORD`, and retags the
wheel (`py3-none-manylinux_2_28_x86_64`, `..._aarch64`, `py3-none-macosx_11_0_universal2`, `py3-none-win_amd64`).
No compiled Python extension is involved: the wheel is the pure package plus one data file that `coal` execs. Linux
binaries are built in the `manylinux_2_28` containers with `libstdc++`/`libgcc` linked statically
(`packaging/wheels/build_linux_binary.sh`), so glibc 2.28 is the only requirement, and `auditwheel show` confirms
the tag; the macOS binary is a universal2 build (`CMAKE_OSX_ARCHITECTURES="arm64;x86_64"`, deployment target 11.0,
checked with `lipo -archs`); the Windows one links the CRT statically (no Visual C++ redistributable needed). To
make one by hand from a local build: `python -m build --wheel` then
`python packaging/wheels/build_platform_wheel.py --binary build/coal-cpp --wheel dist/coal_align-<v>-py3-none-any.whl --platform-tag <tag>`
(`tests/test_wheels.py` does exactly this and installs the result).

## Embedding the library

The algorithm is also a static library, `coal_core` (C++17, standard library only), for programs that build codon
alignments themselves, e.g. a DIAMOND or OrthologR pipeline that already holds the sequences in memory. It has no
global mutable state (the genetic-code tables are a function-local static, built once, thread-safely, and read-only
after that), never calls `exit()`, installs no signal handlers and its public API writes nothing to stdout/stderr: problems are
exceptions (`coal::CoalError` with a stable `code` and the `seq_id`; `coal::InputError` from the file readers).
`thread_one` can be called from several threads at once (coal-cpp does, from its thread pool). Signals, temp files
and output are the `coal-cpp` program's own (`cpp/src/sysio.cpp`, `main.cpp`).

**Installed** (`coal-cpp`, `libcoal_core.a`, the headers under `include/coal/` and a CMake package under
`lib/cmake/coal/`, at the usual GNUInstallDirs places):

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
cmake --install build --prefix /opt/coal
```

```cmake
find_package(coal 0.4 CONFIG REQUIRED)        # -DCMAKE_PREFIX_PATH=/opt/coal; any 0.x >= 0.4 (SameMajorVersion)
target_link_libraries(myapp PRIVATE coal::core)
```

**From source** (`add_subdirectory` or FetchContent): only the library is built, unless you switch on
`COAL_BUILD_CLI` (coal-cpp), `COAL_BUILD_TOOLS` (the test helper) or `COAL_INSTALL` (install rules); as the
top-level project all three are on. The parent's build type and warning flags are left alone.

```cmake
include(FetchContent)
FetchContent_Declare(coal GIT_REPOSITORY https://github.com/plumat-coal/COAL.git GIT_TAG <commit or tag>)
FetchContent_MakeAvailable(coal)              # or: add_subdirectory(third_party/COAL)
target_link_libraries(myapp PRIVATE coal::core)
```

**Headers:** `<coal/core.hpp>` (the algorithm), `<coal/tables.hpp>` (genetic codes), `<coal/seqio.hpp>` (the FASTA /
CLUSTAL readers coal-cpp uses, if you read files) and `<coal/version.hpp>` (`COAL_VERSION_STRING`). The run report
(`report.hpp`, `json.hpp`) is compiled into the library for coal-cpp but is not a public header: its output is tied
to coal-cpp's command line and messages. The API follows the engine and is not frozen before 1.0.

**Minimal API** (all in namespace `coal`; full example: `examples/cpp-embed/`):

```cpp
#include <coal/core.hpp>

const coal::GeneticCode* code = coal::get_code(1);          // NCBI table; nullptr if not supported
std::string row = coal::prepare_row("MKV-AL");              // upper case, '.' -> '-'
std::vector<int> widths(row.size(), 3);                     // one entry per alignment column
coal::update_column_widths(widths, row);                    // for every row: frameshift digits widen a column
coal::Threaded t = coal::thread_one("seq1", row, "ATGAAAGTTGCTCTG", *code);  // throws coal::CoalError
std::string codon_row = coal::expand(row, t.dna, widths);   // "ATGAAAGTT---GCTCTG"
// t.warnings: per-codon mismatch warnings ("pepAlnPos ..."); t.notices: probably-wrong-data notices;
// t.mismatch_cols: columns -nomismatch would remove (coal::keep_slices + coal::mark_nogap_triplets do the removal)
```

Pairing proteins with DNA records (by ID, else by order) is the caller's job; `coal-cpp` does it as PAL2NAL does
(`docs/usage.md`). To check the install: `examples/cpp-embed/` builds against it (`cmake -S examples/cpp-embed -B
build-embed -DCMAKE_PREFIX_PATH=/opt/coal`); CI does this on Linux.
