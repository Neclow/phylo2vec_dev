cmake_minimum_required(VERSION 3.22.1)
project(phylo2vec)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Google Test using FetchContent
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.11.0 # dpkg -l libgtest-dev
)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Find benchmark package
find_package(benchmark REQUIRED)

# Source files
set(SOURCES
    base/to_newick.cpp
    base/to_vector.cpp
    utils/random.cpp
    utils/validation.cpp
    utils/newick.cpp
    utils/vector.cpp
)

set(TEST_SOURCES
    tests/test_main.cpp
    tests/test_v2newick2v.cpp
    tests/test_utils.cpp
)

# Main library
add_library(phylo2vec ${SOURCES})

# Main executable
add_executable(phylo2vec_main main.cpp ${SOURCES})

target_link_libraries(phylo2vec_main phylo2vec)

# Test executable
add_executable(phylo2vec_test ${TEST_SOURCES} ${SOURCES})

# Link against Google Test and Google Mock
target_link_libraries(phylo2vec_test PRIVATE gtest_main)

# Add tests
add_test(NAME phylo2vec_test COMMAND phylo2vec_test)

# Benchmarking
add_executable(phylo2vec_bench benchmarks/bench_base.cpp ${SOURCES})

# Link against Google benchmark
target_link_libraries(phylo2vec_bench PUBLIC benchmark::benchmark)