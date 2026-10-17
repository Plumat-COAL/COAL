cmake_minimum_required(VERSION 3.16)
project(coal_cpp VERSION 0.3.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Threads REQUIRED)

# Core algorithm as a library (no I/O policy), so it can be reused, e.g. in DIAMOND.
add_library(coal_core STATIC cpp/src/tables.cpp cpp/src/core.cpp cpp/src/seqio.cpp)
target_include_directories(coal_core PUBLIC cpp/src)

add_executable(coal-cpp cpp/src/main.cpp)
target_link_libraries(coal-cpp PRIVATE coal_core Threads::Threads)

if(MSVC)
  target_compile_options(coal_core PRIVATE /W4)
  target_compile_options(coal-cpp PRIVATE /W4)
else()
  target_compile_options(coal_core PRIVATE -Wall -Wextra)
  target_compile_options(coal-cpp PRIVATE -Wall -Wextra)
endif()
