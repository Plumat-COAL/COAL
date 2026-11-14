#!/bin/bash
# Bioconda build for coal-align: the C++ engine (coal-cpp) into $PREFIX/bin, then the Python package.
# Linux and macOS only (Bioconda does not build Windows).
set -euxo pipefail

# 1. coal-cpp and coal-cpp-align: C++17, standard library only. COAL_BUILD_TOOLS=OFF skips coal-cpp-dump (a test
#    helper). Only the two programs are installed: this package is the command line and the Python package. The
#    static library libcoal_core.a, its headers and CMake package (what `cmake --install` adds for C++ projects
#    that embed the engine) are not shipped: nothing in Bioconda links against them, and a static library in a
#    conda package needs a reason. COAL_INSTALL=OFF: no install rules are generated.
#    macOS: Bioconda's osx-64 deployment target (10.13) predates std::filesystem in Apple's libc++; conda ships its
#    own modern libc++, so the availability checks are switched off (conda-forge knowledge base, "Newer C++
#    features with old SDK").
if [[ "${target_platform:-}" == osx-* ]]; then
  export CXXFLAGS="${CXXFLAGS:-} -D_LIBCPP_DISABLE_AVAILABILITY"
fi
cmake -S . -B build-conda \
  -DCMAKE_BUILD_TYPE=Release \
  -DCOAL_BUILD_TOOLS=OFF \
  -DCOAL_INSTALL=OFF \
  ${CMAKE_ARGS:-}
cmake --build build-conda --config Release -j "${CPU_COUNT:-2}"

mkdir -p "${PREFIX}/bin"
install -m 0755 build-conda/coal-cpp "${PREFIX}/bin/coal-cpp"
install -m 0755 build-conda/coal-cpp-align "${PREFIX}/bin/coal-cpp-align"
rm -rf build-conda

# 2. The Python engine: `coal`, `pal2nal`, `import coal`; examples under share/coal/examples. No dependencies.
"${PYTHON}" -m pip install . --no-deps --no-build-isolation --no-cache-dir -vv

# `coal` hands inputs >= 1 MB to coal-cpp when it finds one of its own version next to the Python executable
# ($PREFIX/bin, where both land here) or on PATH, so no COAL_CPP setting is needed in a conda environment.
# (`coal --version` prints "coal X.Y.Z", `coal-cpp --version` prints "coal X.Y.Z (C++)".)
test "$("${PREFIX}/bin/coal" --version | cut -d' ' -f2)" = "$("${PREFIX}/bin/coal-cpp" --version | cut -d' ' -f2)"
