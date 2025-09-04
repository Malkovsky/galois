cmake_minimum_required(VERSION 3.12)
project(gf256)

set(CMAKE_CXX_FLAGS "-march=native")
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.9.4  # Replace with the desired version/tag
)
FetchContent_MakeAvailable(googlebenchmark)
add_executable(benchmarks
    src/field.cc
    third_party/gf256/gf256.cpp
    benchmarks/matrix_multiplication.cc)
target_include_directories(benchmarks
    PUBLIC src
    PUBLIC third_party
)
set_property(TARGET benchmarks PROPERTY CXX_STANDARD 20)

target_link_libraries(benchmarks
    benchmark
    benchmark_main)
add_test(
    NAME Benchmarks
    COMMAND benchmarks --benchmark_out=bm_report.csv --benchmark_out_format=csv
)

add_executable(gf_unittests
    src/field.cc
    third_party/gf256/gf256.cpp
    tests/field_tests.cc)
target_include_directories(gf_unittests
    PUBLIC ${GOOGLETEST_SOURCE_DIR}/src
    PUBLIC src
    PUBLIC third_party
)
set_property(TARGET gf_unittests PROPERTY CXX_STANDARD 20)
target_link_libraries(gf_unittests
    gtest
    gtest_main)
add_test(
    NAME Unittests
    COMMAND gf_unittests
)
enable_testing()