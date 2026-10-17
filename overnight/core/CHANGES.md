# C++ core: changes (branch overnight/core)

One row per change. Reference behaviour: the reviewed Python engine on `origin/keanu_overnight/integration`
(read with `git show`, not merged). IDs appear at the start of the commit messages.

| ID | Where | Before | After | Why | Test that proves it |
|---|---|---|---|---|---|
| CPP-1 | CMakeLists.txt | `/W4` (MSVC) or `-Wall -Wextra`; no way to fail the build on warnings | Options `COAL_WARNINGS_AS_ERRORS` (adds `/WX` or `-Werror`) and `COAL_BUILD_TOOLS` (test helper `coal-cpp-dump`); MSVC `/W4 /permissive- /utf-8 /EHsc` + `_CRT_SECURE_NO_WARNINGS NOMINMAX`; GCC/Clang `-Wall -Wextra -Wpedantic -Wshadow` | Brief item 1: CI-friendly warnings-as-errors on Windows, macOS and Linux | `cmake -DCOAL_WARNINGS_AS_ERRORS=ON` builds with GCC 11 (Linux); clang 19 `-Werror` compiles every file for x86_64-windows-gnu, aarch64-macos and x86_64-linux-gnu. MSVC itself not available here (see LOG) |
