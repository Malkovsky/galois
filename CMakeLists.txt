cmake_minimum_required(VERSION 3.12)
project(gf256)

#
# Compiler options
#
# option(${PROJECT_NAME}_WARNINGS_AS_ERRORS "Treat compiler warnings as errors." OFF)
# option(${PROJECT_NAME}_USE_MSVC_MultiThreadDLL "Link dynamically-link mutithread library. This flag only take effect when it is MSVC" OFF)

# if(MSVC)
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMTD")
# endif()

# # Add this to your CMakeLists.txt
# if(MSVC)
#   # Use dynamic runtime libraries consistently
#   set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDebugDLL")
  
#   # For release builds:
#   # set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
# endif()

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

add_executable(gf_unittests
    src/field.cc
    tests/field_tests.cc)
target_include_directories(gf_unittests
    PUBLIC ${GOOGLETEST_SOURCE_DIR}/src
    PUBLIC src
)
set_property(TARGET gf_unittests PROPERTY CXX_STANDARD 20)
target_link_libraries(gf_unittests
    gtest
    gtest_main)
add_test(
    NAME Unittests
    COMMAND gf_unittests
)