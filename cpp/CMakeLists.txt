cmake_minimum_required(VERSION 3.22.1)
project(
    phylo2vec
    VERSION 1.0.0
    LANGUAGES CXX
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include(FetchContent)

# Include fmt using FetchContent
# FetchContent_Declare(
#   fmt
#   GIT_REPOSITORY https://github.com/fmtlib/fmt
#   GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281) # 10.2.1
# FetchContent_MakeAvailable(fmt)

# Source files
set(SOURCES
    base/avl.cpp
    base/to_newick.cpp
    base/to_vector.cpp
    matrix/to_newick.cpp
    utils/random.cpp
    utils/validation.cpp
    utils/newick.cpp
    utils/vector.cpp
)

set(BENCH_SOURCES
    benchmarks/bench_base.cpp
)

set(TEST_SOURCES
    tests/test_main.cpp
    tests/test_v2newick2v.cpp
    tests/test_utils.cpp
)

# Main library
add_library(phylo2vec_cpp STATIC ${SOURCES})

# Add git submodules
add_subdirectory(extern)

# Main executable
# add_executable(phylo2vec_main main.cpp ${SOURCES})

# target_link_libraries(phylo2vec_main phylo2vec_cpp)

if(BUILD_TESTING)
    # Include Google Test using FetchContent
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.11.0 # dpkg -l libgtest-dev
    )

    FetchContent_MakeAvailable(googletest)
    # Test executable
    add_executable(phylo2vec_test ${TEST_SOURCES} ${SOURCES})

    # Link against Google Test and Google Mock
    target_link_libraries(phylo2vec_test PRIVATE gtest_main)

    # Add tests
    add_test(NAME phylo2vec_test COMMAND phylo2vec_test)
endif()

# Benchmarking
if(BUILD_BENCHMARKS)
    # Find benchmark package
    find_package(benchmark REQUIRED)

    # Add bench_base.cpp
    add_executable(phylo2vec_bench ${BENCH_SOURCES} ${SOURCES})

    # Link against Google benchmark
    target_link_libraries(phylo2vec_bench PUBLIC benchmark::benchmark)
endif()

if(BUILD_PYTHON)
  add_subdirectory(python)
endif()