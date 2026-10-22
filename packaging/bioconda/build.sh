#!/bin/bash
# Bioconda build for coal-align: the C++ engine (coal-cpp) into $PREFIX/bin, then the Python package.
# Linux and macOS only (Bioconda does not build Windows).
set -euxo pipefail

# 1. coal-cpp: C++17, standard library only. COAL_BUILD_TOOLS=OFF skips coal-cpp-dump (a test helper).
cmake -S . -B build-conda \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DCOAL_BUILD_TOOLS=OFF \
  ${CMAKE_ARGS:-}
cmake --build build-conda --config Release -j "${CPU_COUNT:-2}"

mkdir -p "${PREFIX}/bin"
if grep -q 'install(TARGETS coal-cpp' CMakeLists.txt; then
  # CMakeLists.txt has install rules (COAL_INSTALL, from the feat/stdin-and-cpp-lib branch): installs bin/coal-cpp and,
  # for C++ users embedding the engine, lib/libcoal_core.a, include/coal/ and lib/cmake/coal/.
  cmake --install build-conda --config Release
else
  # Older CMakeLists.txt without install rules (0.4.2 and earlier): copy the binary by hand.
  install -m 0755 build-conda/coal-cpp "${PREFIX}/bin/coal-cpp"
fi
rm -rf build-conda

# 2. The Python engine: `coal`, `pal2nal`, `import coal`; examples under share/coal/examples. No dependencies.
"${PYTHON}" -m pip install . --no-deps --no-build-isolation --no-cache-dir -vv

# `coal` hands inputs >= 1 MB to coal-cpp when it finds one of its own version next to the Python executable
# ($PREFIX/bin, where both land here) or on PATH, so no COAL_CPP setting is needed in a conda environment.
# (`coal --version` prints "coal X.Y.Z", `coal-cpp --version` prints "coal X.Y.Z (C++)".)
test "$("${PREFIX}/bin/coal" --version | cut -d' ' -f2)" = "$("${PREFIX}/bin/coal-cpp" --version | cut -d' ' -f2)"
