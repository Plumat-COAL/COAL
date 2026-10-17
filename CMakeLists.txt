cmake_minimum_required(VERSION 3.16)
project(coal_cpp VERSION 0.4.0 LANGUAGES CXX)

# coal-cpp: C++17, standard library only. Builds with MSVC, Clang and GCC (Windows, macOS, Linux).
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
# CI: add -DCOAL_WARNINGS_AS_ERRORS=ON (warnings fail the build).

option(COAL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (for CI)" OFF)
option(COAL_BUILD_TOOLS "Build coal-cpp-dump (test helper used by tests/test_cpp_*.py)" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Threads REQUIRED)

# Core algorithm as a library (no I/O policy), so it can be reused, e.g. in DIAMOND.
add_library(coal_core STATIC cpp/src/tables.cpp cpp/src/core.cpp cpp/src/seqio.cpp cpp/src/json.cpp cpp/src/report.cpp)
target_include_directories(coal_core PUBLIC cpp/src)

add_executable(coal-cpp cpp/src/main.cpp cpp/src/args.cpp)
target_link_libraries(coal-cpp PRIVATE coal_core Threads::Threads)

set(COAL_TARGETS coal_core coal-cpp)
if(COAL_BUILD_TOOLS)
  add_executable(coal-cpp-dump cpp/tools/dump.cpp)
  target_link_libraries(coal-cpp-dump PRIVATE coal_core Threads::Threads)
  list(APPEND COAL_TARGETS coal-cpp-dump)
endif()

foreach(t IN LISTS COAL_TARGETS)
  if(MSVC)
    # /utf-8: sources are UTF-8; _CRT_SECURE_NO_WARNINGS: std::getenv/strerror are used on purpose
    target_compile_options(${t} PRIVATE /W4 /permissive- /utf-8 /EHsc)
    target_compile_definitions(${t} PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
    if(COAL_WARNINGS_AS_ERRORS)
      target_compile_options(${t} PRIVATE /WX)
    endif()
  else()
    target_compile_options(${t} PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
    if(COAL_WARNINGS_AS_ERRORS)
      target_compile_options(${t} PRIVATE -Werror)
    endif()
  endif()
endforeach()
