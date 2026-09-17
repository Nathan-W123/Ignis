# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.16)
project(Ignis VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(IGNIS_BUILD_TESTS  "Build the Catch2 test suite"     ON)
option(IGNIS_BUILD_APPS   "Build the command-line tools"    ON)
option(IGNIS_BUILD_BENCH  "Build the benchmark driver"      ON)
option(IGNIS_WERROR       "Treat compiler warnings as errors" OFF)
option(IGNIS_ENABLE_LTO   "Enable interprocedural optimisation in Release" OFF)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
include(FetchContent)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found locally -- fetching release 3.4.0")
  FetchContent_Declare(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(Eigen3)
endif()

find_package(yaml-cpp QUIET)
if(NOT yaml-cpp_FOUND AND NOT TARGET yaml-cpp::yaml-cpp)
  message(STATUS "yaml-cpp not found locally -- fetching 0.8.0")
  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        0.8.0
    GIT_SHALLOW    TRUE)
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
  set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(yaml-cpp)
endif()
if(TARGET yaml-cpp::yaml-cpp)
  set(IGNIS_YAML_TARGET yaml-cpp::yaml-cpp)
else()
  set(IGNIS_YAML_TARGET yaml-cpp)
endif()

find_package(Threads REQUIRED)

# ---------------------------------------------------------------------------
# Version header
# ---------------------------------------------------------------------------
find_package(Git QUIET)
set(IGNIS_GIT_DESCRIBE "unknown")
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE IGNIS_GIT_DESCRIBE
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
endif()
set(IGNIS_DATA_DIR "${PROJECT_SOURCE_DIR}/data" CACHE PATH "Directory holding Ignis data files")
configure_file(${PROJECT_SOURCE_DIR}/include/ignis/core/Version.hpp.in
               ${PROJECT_BINARY_DIR}/generated/ignis/core/Version.hpp @ONLY)

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
file(GLOB_RECURSE IGNIS_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(ignis STATIC ${IGNIS_SOURCES})
add_library(Ignis::ignis ALIAS ignis)

target_include_directories(ignis PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(ignis PUBLIC Eigen3::Eigen Threads::Threads PRIVATE ${IGNIS_YAML_TARGET})

target_compile_options(ignis PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
    -Wall -Wextra -Wpedantic -Wshadow -Wno-unused-parameter>
  $<$<BOOL:${IGNIS_WERROR}>:-Werror>)

if(IGNIS_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
  if(ipo_ok)
    set_property(TARGET ignis PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  endif()
endif()

# ---------------------------------------------------------------------------
# Applications
# ---------------------------------------------------------------------------
if(IGNIS_BUILD_APPS)
  foreach(app ignis_equilibrium ignis_engine ignis_nozzle ignis_transient ignis_sweep ignis_mc)
    add_executable(${app} ${PROJECT_SOURCE_DIR}/apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE ignis ${IGNIS_YAML_TARGET})
    set_target_properties(${app} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
  endforeach()
endif()

if(IGNIS_BUILD_BENCH)
  add_executable(ignis_bench ${PROJECT_SOURCE_DIR}/bench/ignis_bench.cpp)
  target_link_libraries(ignis_bench PRIVATE ignis ${IGNIS_YAML_TARGET})
  set_target_properties(ignis_bench PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(IGNIS_BUILD_TESTS)
  enable_testing()
  find_package(Catch2 3 QUIET)
  if(NOT Catch2_FOUND)
    message(STATUS "Catch2 v3 not found locally -- fetching v3.5.2")
    FetchContent_Declare(Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v3.5.2
      GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  else()
    list(APPEND CMAKE_MODULE_PATH ${Catch2_DIR})
  endif()
  add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "Ignis ${PROJECT_VERSION}  (${IGNIS_GIT_DESCRIBE})")
message(STATUS "  build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "  compiler        : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  data directory  : ${IGNIS_DATA_DIR}")
message(STATUS "  tests / apps    : ${IGNIS_BUILD_TESTS} / ${IGNIS_BUILD_APPS}")
message(STATUS "")
