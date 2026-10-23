#!/usr/bin/env bash
# Build coal-cpp for a manylinux_2_28 wheel. Run inside quay.io/pypa/manylinux_2_28_x86_64 or
# quay.io/pypa/manylinux_2_28_aarch64 (release.yml does: docker run -v "$PWD:/io" <image> bash /io/packaging/wheels/build_linux_binary.sh),
# whose glibc is 2.28 and whose GCC 14 toolset links libstdc++ and libgcc statically here, so the binary needs nothing
# but glibc >= 2.28 (RHEL 8 / Debian 10 / Ubuntu 18.10 and later). auditwheel show on the wheel confirms the tag.
#
#   build_linux_binary.sh [SOURCE_DIR=/io] [BUILD_DIR=SOURCE_DIR/build-manylinux]
set -euo pipefail

src=${1:-/io}
out=${2:-$src/build-manylinux}

if ! command -v cmake > /dev/null 2>&1; then  # the images ship cmake via pipx; a fallback for a plain el8 container
    pipx install cmake > /dev/null 2>&1 || /opt/python/cp312-cp312/bin/python -m pip install --quiet cmake
    export PATH="$HOME/.local/bin:/opt/python/cp312-cp312/bin:$PATH"
fi

cmake -S "$src" -B "$out" -DCMAKE_BUILD_TYPE=Release -DCOAL_WARNINGS_AS_ERRORS=ON -DCOAL_BUILD_TOOLS=OFF \
    -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc"
cmake --build "$out" -j "$(nproc)"

bin="$out/coal-cpp"
test -x "$bin"
"$bin" --version
ldd "$bin"
# Nothing but glibc may be dynamic (libstdc++ and libgcc_s are linked in); auditwheel would otherwise demand a newer tag.
extra=$(ldd "$bin" | grep -vE 'linux-vdso|ld-linux|ld64|libc\.so|libm\.so|libpthread\.so|libdl\.so|librt\.so|^[[:space:]]*$' || true)
if [ -n "$extra" ]; then
    echo "coal-cpp links shared libraries other than glibc:" >&2
    echo "$extra" >&2
    exit 1
fi
echo "built $bin ($(uname -m)), glibc only"
