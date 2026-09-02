cmake_minimum_required(VERSION 3.24)
project(gf256 LANGUAGES C CXX)

execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GF256_GIT_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT GF256_GIT_REVISION)
  set(GF256_GIT_REVISION unknown)
endif()
execute_process(
  COMMAND git status --porcelain
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GF256_GIT_DIRTY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(GF256_GIT_DIRTY)
  string(APPEND GF256_GIT_REVISION "-dirty")
endif()

option(GF256_BUILD_TESTS "Build unit tests" ON)
option(GF256_BUILD_BENCHMARKS "Build benchmarks" ON)
option(GF256_BUILD_REFERENCE_BENCHMARKS "Build native RS baselines" OFF)
option(GF256_ENABLE_NATIVE_ISA "Compile owned sources for the build host" OFF)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.9.4
)

add_library(gf256_core STATIC
    src/field.cc
    src/matrix.cc
    src/reed_solomon/common.cc
    src/reed_solomon/leopard.cc
    src/reed_solomon/xdrs.cc
    src/reed_solomon/lin_chung_han/transform.cc
    src/reed_solomon/lin_chung_han/kernels.cc
    src/reed_solomon/lin_chung_han/gfni.cc
)
target_compile_features(gf256_core PUBLIC cxx_std_20)
target_include_directories(gf256_core
  PUBLIC ${PROJECT_SOURCE_DIR}/include
  PRIVATE ${PROJECT_SOURCE_DIR}/src)

if(GF256_ENABLE_NATIVE_ISA)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gf256_core PRIVATE -march=native)
  else()
    message(WARNING "GF256_ENABLE_NATIVE_ISA has no flags for this compiler")
  endif()
endif()

if(GF256_BUILD_TESTS)
  FetchContent_MakeAvailable(googletest)
  enable_testing()

  add_executable(gf_unittests tests/field_tests.cc)
  target_link_libraries(gf_unittests PRIVATE gf256_core gtest gtest_main)
  add_test(NAME Unittests COMMAND gf_unittests)

  add_executable(lch_rs_unittests tests/lch_tests.cc tests/rs_tests.cc)
  target_include_directories(lch_rs_unittests PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(lch_rs_unittests PRIVATE gf256_core gtest gtest_main)
  add_test(NAME LCHRSUnittests COMMAND lch_rs_unittests)
endif()

if(GF256_BUILD_BENCHMARKS)
  FetchContent_MakeAvailable(googlebenchmark)

  add_executable(benchmarks benchmarks/matrix_multiplication.cc)
  target_link_libraries(benchmarks PRIVATE gf256_core benchmark benchmark_main)

  add_executable(rs_verbose_benchmarks benchmarks/lch_rs_benchmarks.cc)
  add_executable(rs_benchmarks benchmarks/lch_rs_benchmarks.cc)
  target_compile_definitions(rs_benchmarks PRIVATE
    GF256_HIGH_LEVEL_BENCHMARKS_ONLY=1)

  foreach(benchmark_target IN ITEMS rs_verbose_benchmarks rs_benchmarks)
    target_link_libraries(${benchmark_target} PRIVATE gf256_core benchmark)
    target_compile_definitions(
      ${benchmark_target} PRIVATE
        GF256_GIT_REVISION="${GF256_GIT_REVISION}"
        GF256_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
        GF256_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
        GF256_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}")
  endforeach()

  if(GF256_ENABLE_NATIVE_ISA AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(benchmark_target IN ITEMS rs_verbose_benchmarks rs_benchmarks)
      target_compile_options(${benchmark_target} PRIVATE -march=native)
      target_compile_definitions(${benchmark_target} PRIVATE
        GF256_NATIVE_ISA="-march=native")
    endforeach()
  else()
    foreach(benchmark_target IN ITEMS rs_verbose_benchmarks rs_benchmarks)
      target_compile_definitions(${benchmark_target} PRIVATE
        GF256_NATIVE_ISA="disabled")
    endforeach()
  endif()

  if(GF256_BUILD_REFERENCE_BENCHMARKS)
    add_library(gf256_leopard_reference STATIC
      third_party/leopard/leopard.cpp
      third_party/leopard/LeopardCommon.cpp
      third_party/leopard/LeopardFF8.cpp
      third_party/leopard/LeopardFF16.cpp)
    target_include_directories(gf256_leopard_reference PUBLIC
      third_party/leopard)
    target_compile_features(gf256_leopard_reference PRIVATE cxx_std_11)
    target_compile_options(gf256_leopard_reference PRIVATE -march=native)

    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" OR
       NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64|i.86")
      message(FATAL_ERROR
        "GF256_BUILD_REFERENCE_BENCHMARKS requires an x86 GNU/Clang build for XDRS")
    endif()
    add_library(gf256_xdrs_reference STATIC
      third_party/xdrs/src/P_function.cpp)
    target_include_directories(gf256_xdrs_reference PUBLIC
      third_party/xdrs/src)
    target_compile_options(gf256_xdrs_reference PRIVATE
      -mavx2 -include cstring)

    foreach(benchmark_target IN ITEMS rs_verbose_benchmarks rs_benchmarks)
      target_sources(${benchmark_target} PRIVATE
        benchmarks/reference/leopard_benchmarks.cc
        benchmarks/reference/xdrs_benchmarks.cc)
      target_link_libraries(${benchmark_target} PRIVATE
        gf256_leopard_reference gf256_xdrs_reference)
      target_compile_definitions(${benchmark_target} PRIVATE
        GF256_HAVE_REFERENCE_BENCHMARKS=1)
    endforeach()
  endif()
endif()
