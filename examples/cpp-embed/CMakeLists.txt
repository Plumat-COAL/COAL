# Example: a program that embeds the COAL library (coal::core) to codon-align sequences held in memory.
#
#   cmake -S . -B build -DCMAKE_PREFIX_PATH=/where/coal/is/installed
#   cmake --build build && ./build/coal-embed-example
#
# COAL installed with:  cmake -S <COAL> -B <COAL>/build && cmake --install <COAL>/build --prefix /where/...
cmake_minimum_required(VERSION 3.16)
project(coal_embed_example LANGUAGES CXX)

find_package(coal 0.4 CONFIG REQUIRED)

add_executable(coal-embed-example main.cpp)
target_link_libraries(coal-embed-example PRIVATE coal::core)
