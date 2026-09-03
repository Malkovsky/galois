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
option(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT
  "Build retained GFNI512 radix-8 tests and verbose benchmarks" OFF)

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
    src/lin_chung_han/transform.cc
    src/lin_chung_han/kernels.cc
    src/lin_chung_han/gfni.cc
    src/reed_solomon/code_parameters.cc
    src/reed_solomon/decoder/shared.cc
    src/reed_solomon/decoder/low_rate.cc
    src/reed_solomon/decoder/high_rate.cc
    src/reed_solomon/lch_encoder.cc
    src/reed_solomon/lch_decoder.cc
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

if(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
  if(NOT GF256_ENABLE_NATIVE_ISA)
    message(FATAL_ERROR
      "GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT requires GF256_ENABLE_NATIVE_ISA")
  endif()
  add_library(gf256_gfni512_radix8_experiment STATIC
    src/lin_chung_han/experiment/gfni512_radix8.cc
    src/reed_solomon/experiment/gfni512_radix8.cc)
  target_compile_features(gf256_gfni512_radix8_experiment PRIVATE cxx_std_20)
  target_include_directories(gf256_gfni512_radix8_experiment PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(gf256_gfni512_radix8_experiment PUBLIC gf256_core)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gf256_gfni512_radix8_experiment PRIVATE
      -march=native)
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
  if(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
    target_link_libraries(lch_rs_unittests PRIVATE
      gf256_gfni512_radix8_experiment)
    target_compile_definitions(lch_rs_unittests PRIVATE
      GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT=1)
  endif()
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
  if(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
    target_compile_definitions(rs_verbose_benchmarks PRIVATE
      GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT=1)
    target_link_libraries(rs_verbose_benchmarks PRIVATE
      gf256_gfni512_radix8_experiment)
  endif()

  foreach(benchmark_target IN ITEMS rs_verbose_benchmarks rs_benchmarks)
    target_include_directories(${benchmark_target} PRIVATE
      ${PROJECT_SOURCE_DIR}/src)
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
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" OR
       NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
      message(FATAL_ERROR
        "GF256_BUILD_REFERENCE_BENCHMARKS requires an x86 GNU/Clang build")
    endif()
    find_program(GF256_NASM_EXECUTABLE nasm)
    if(NOT GF256_NASM_EXECUTABLE)
      message(FATAL_ERROR
        "GF256_BUILD_REFERENCE_BENCHMARKS requires NASM for ISA-L")
    endif()
    set(ISAL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ISAL_BUILD_PERF_TESTS OFF CACHE BOOL "" FORCE)
    set(ISAL_BUILD_FUZZ_TESTS OFF CACHE BOOL "" FORCE)
    set(ISAL_BUILD_ISAL_SHIM OFF CACHE BOOL "" FORCE)
    set(ISAL_BUILD_IGZIP_CLI OFF CACHE BOOL "" FORCE)
    get_property(GF256_BUILD_SHARED_LIBS_CACHE_TYPE CACHE BUILD_SHARED_LIBS
      PROPERTY TYPE)
    if(GF256_BUILD_SHARED_LIBS_CACHE_TYPE)
      set(GF256_SAVED_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    endif()
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    add_subdirectory(third_party/isa-l third_party/isa-l-build
      EXCLUDE_FROM_ALL)
    if(GF256_BUILD_SHARED_LIBS_CACHE_TYPE)
      set(BUILD_SHARED_LIBS ${GF256_SAVED_BUILD_SHARED_LIBS}
        CACHE BOOL "" FORCE)
    else()
      unset(BUILD_SHARED_LIBS CACHE)
    endif()

    add_library(gf256_leopard_reference STATIC
      third_party/leopard/leopard.cpp
      third_party/leopard/LeopardCommon.cpp
      third_party/leopard/LeopardFF8.cpp
      third_party/leopard/LeopardFF16.cpp)
    target_include_directories(gf256_leopard_reference PUBLIC
      third_party/leopard)
    target_compile_features(gf256_leopard_reference PRIVATE cxx_std_11)
    target_compile_options(gf256_leopard_reference PRIVATE -march=native)

    add_library(gf256_xdrs_reference STATIC
      third_party/xdrs/src/P_function.cpp)
    target_include_directories(gf256_xdrs_reference PUBLIC
      third_party/xdrs/src)
    target_compile_options(gf256_xdrs_reference PRIVATE
      -mavx2 -include cstring)

    foreach(benchmark_target IN ITEMS rs_verbose_benchmarks rs_benchmarks)
      target_sources(${benchmark_target} PRIVATE
        benchmarks/reference/isal_benchmarks.cc
        benchmarks/reference/leopard_benchmarks.cc
        benchmarks/reference/xdrs_benchmarks.cc)
      target_link_libraries(${benchmark_target} PRIVATE
        gf256_leopard_reference gf256_xdrs_reference isal)
      target_compile_definitions(${benchmark_target} PRIVATE
        GF256_HAVE_REFERENCE_BENCHMARKS=1)
    endforeach()
  endif()
endif()
