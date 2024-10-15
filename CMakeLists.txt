#
# Copyright (C) 2009-2022 The ESPResSo project
# Copyright (C) 2009,2010
#   Max-Planck-Institute for Polymer Research, Theory Group
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.25.1)
cmake_policy(VERSION 3.25.1)
if(POLICY CMP0167)
  # use BoostConfig.cmake shipped with Boost 1.70+ instead of the one in CMake
  cmake_policy(SET CMP0167 NEW)
endif()
message(STATUS "CMake version: ${CMAKE_VERSION}")
# CMake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# project info
project(ESPResSo)
set(PROJECT_VERSION "4.3-dev")

# C++ standard
enable_language(CXX)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# require minimal compiler versions
macro(espresso_minimal_compiler_version)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "${ARGV0}" AND CMAKE_CXX_COMPILER_VERSION
                                                   VERSION_LESS "${ARGV1}")
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID} \
                 ${CMAKE_CXX_COMPILER_VERSION} (required version >= ${ARGV1})")
  endif()
endmacro()

espresso_minimal_compiler_version("GNU" 10.5.0)
espresso_minimal_compiler_version("Clang" 14.0.0)
espresso_minimal_compiler_version("AppleClang" 14.0.0)
espresso_minimal_compiler_version("Intel" 2021.9)
espresso_minimal_compiler_version("IntelLLVM" 2023.1)

include(FeatureSummary)
include(GNUInstallDirs)
include(FetchContent)
include(espresso_option_enum)
include(espresso_enable_avx2_support)

#
# CMake internal vars
#

# Select the build type
espresso_option_enum(
  varname "CMAKE_BUILD_TYPE" help_text "build type" default_value "Release"
  possible_values
  "Debug;Release;RelWithDebInfo;MinSizeRel;Coverage;RelWithAssert")
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE} -Og -g")
set(CMAKE_CXX_FLAGS_RELWITHASSERT "${CMAKE_CXX_FLAGS_RELWITHASSERT} -O3 -g")

# On Mac OS X, first look for other packages, then frameworks
set(CMAKE_FIND_FRAMEWORK LAST)

# avoid applying patches twice with FetchContent
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

#
# External dependencies declarations
#

# cmake-format: off
FetchContent_Declare(
  walberla
  GIT_REPOSITORY https://i10git.cs.fau.de/walberla/walberla.git
  GIT_TAG        b0842e1a493ce19ef1bbb8d2cf382fc343970a7f
)
FetchContent_Declare(
  stokesian_dynamics
  GIT_REPOSITORY https://github.com/hmenke/espresso-stokesian-dynamics.git
  GIT_TAG        862a7537a366f0c32f0c25e46bd107bea590faea
)
FetchContent_Declare(
  caliper
  GIT_REPOSITORY https://github.com/LLNL/Caliper.git
  GIT_TAG        v2.10.0
)
# cmake-format: on

# ##############################################################################
# User input options
# ##############################################################################

option(ESPRESSO_BUILD_WITH_PYTHON "Build with Python bindings" ON)
option(ESPRESSO_BUILD_WITH_GSL "Build with GSL support" OFF)
option(ESPRESSO_BUILD_WITH_FFTW "Build with FFTW support" ON)
option(ESPRESSO_BUILD_WITH_CUDA "Build with GPU support" OFF)
option(ESPRESSO_BUILD_WITH_HDF5 "Build with HDF5 support" OFF)
option(ESPRESSO_BUILD_TESTS "Enable tests" ON)
option(ESPRESSO_BUILD_WITH_SCAFACOS "Build with ScaFaCoS support" OFF)
option(ESPRESSO_BUILD_WITH_STOKESIAN_DYNAMICS "Build with Stokesian Dynamics"
       OFF)
option(ESPRESSO_BUILD_WITH_WALBERLA
       "Build with waLBerla lattice-Boltzmann support" OFF)
option(ESPRESSO_BUILD_WITH_WALBERLA_AVX
       "Build waLBerla lattice-Boltzmann with AVX vectorization" OFF)
option(ESPRESSO_BUILD_WITH_WALBERLA_FFT "Build waLBerla with FFT support" OFF)
option(ESPRESSO_BUILD_BENCHMARKS "Enable benchmarks" OFF)
option(ESPRESSO_BUILD_WITH_VALGRIND "Build with Valgrind instrumentation" OFF)
option(ESPRESSO_BUILD_WITH_CALIPER "Build with Caliper instrumentation" OFF)
option(ESPRESSO_BUILD_WITH_CPPCHECK "Run Cppcheck during compilation" OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  option(ESPRESSO_BUILD_WITH_CLANG_TIDY "Run Clang-Tidy during compilation" OFF)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                            "GNU")
  option(ESPRESSO_BUILD_WITH_COVERAGE
         "Generate code coverage report for C++ code" OFF)
  option(ESPRESSO_BUILD_WITH_COVERAGE_PYTHON
         "Generate code coverage report for Python code" OFF)
  option(ESPRESSO_BUILD_WITH_ASAN "Build with address sanitizer" OFF)
  option(ESPRESSO_BUILD_WITH_UBSAN "Build with undefined behavior sanitizer"
         OFF)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT APPLE)
  option(
    ESPRESSO_BUILD_WITH_MSAN
    "Build with memory sanitizer (experimental; requires a memory-sanitized Python interpreter)"
    OFF)
endif()
option(
  ESPRESSO_ADD_OMPI_SINGLETON_WARNING
  "Add a runtime warning in the pypresso script for NUMA architectures that aren't supported in singleton mode by Open MPI 4.x"
  ON)
option(ESPRESSO_WARNINGS_ARE_ERRORS
       "Treat warnings as errors during compilation" OFF)
option(ESPRESSO_BUILD_WITH_CCACHE "Use ccache compiler invocation." OFF)
option(ESPRESSO_INSIDE_DOCKER "Set this to ON when running inside Docker." OFF)
mark_as_advanced(ESPRESSO_INSIDE_DOCKER)
set(ESPRESSO_TEST_TIMEOUT "300"
    CACHE STRING "Timeout in seconds for each testsuite test")

if(ESPRESSO_BUILD_WITH_CCACHE)
  find_program(CCACHE_PROGRAM ccache REQUIRED)
  if(CCACHE_PROGRAM)
    message(STATUS "Found ccache: ${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_CUDA_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
  endif()
endif()

# Write compile commands to file, for various tools...
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# choose the name of the config file
set(ESPRESSO_MYCONFIG_NAME "myconfig.hpp"
    CACHE STRING "Default name of the local config file")

# Check which config file to use
include(espresso_myconfig)

#
# Pretty function
#

include(CheckCXXSourceCompiles)

# cross-platform macro to print the function name in error messages
set(ESPRESSO_PRETTY_FUNCTION_EXTENSION __func__)

# search for a supported compiler extension that prints the function name as
# well as its list of arguments, return type and namespace
foreach(func_name __PRETTY_FUNCTION__ __FUNCSIG__ __FUNCTION__)
  check_cxx_source_compiles(
    "
     #include <string>
     int main() { std::string(${func_name}); }
     " result${func_name})
  if(result${func_name})
    set(ESPRESSO_PRETTY_FUNCTION_EXTENSION ${func_name})
    break()
  endif(result${func_name})
endforeach()

#
# AVX2 support
#

include(CheckCXXCompilerFlag)

add_library(espresso_avx_flags INTERFACE)
add_library(espresso::avx_flags ALIAS espresso_avx_flags)

#
# Interface libraries
#

# CUDA compiler
if(ESPRESSO_BUILD_WITH_CUDA)
  cmake_minimum_required(VERSION 3.25.2)
  include(CheckLanguage)
  enable_language(CUDA)
  check_language(CUDA)
  espresso_minimal_compiler_version("GNU" 11.4.0)
  espresso_minimal_compiler_version("Clang" 17.0.0)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(ESPRESSO_MINIMAL_CUDA_VERSION 12.0)
  find_package(CUDAToolkit ${ESPRESSO_MINIMAL_CUDA_VERSION} REQUIRED)
  if(NOT DEFINED ESPRESSO_CMAKE_CUDA_ARCHITECTURES)
    if("$ENV{CUDAARCHS}" STREQUAL "")
      # 1. sm_61: GTX-1000 series (Pascal)
      # 2. sm_75: RTX-2000 series (Turing)
      # 3. sm_86: RTX-3000 series (Ampere)
      # 4. sm_89: RTX-4000 series (Ada)
      set(ESPRESSO_CUDA_ARCHITECTURES "61;75")
    else()
      set(ESPRESSO_CUDA_ARCHITECTURES "$ENV{CUDAARCHS}")
    endif()
    set(ESPRESSO_CMAKE_CUDA_ARCHITECTURES "${ESPRESSO_CUDA_ARCHITECTURES}"
        CACHE INTERNAL "")
  endif()
  set(CMAKE_CUDA_ARCHITECTURES "${ESPRESSO_CMAKE_CUDA_ARCHITECTURES}")
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    find_package(CUDACompilerNVCC ${ESPRESSO_MINIMAL_CUDA_VERSION} REQUIRED)
  elseif(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
    if(ESPRESSO_BUILD_WITH_COVERAGE)
      message(
        FATAL_ERROR
          "Cannot enable code coverage with Clang as the CUDA compiler")
    endif()
    find_package(CUDACompilerClang 17.0 REQUIRED)
  else()
    message(FATAL_ERROR "Unknown CUDA compiler '${CMAKE_CUDA_COMPILER_ID}'")
  endif()
endif()

# Python interpreter and Cython interface library
if(ESPRESSO_BUILD_WITH_PYTHON)
  find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Cython 0.29.28...<3.0.10 REQUIRED)
  find_program(IPYTHON_EXECUTABLE NAMES jupyter ipython3 ipython)
endif()

#
# Installation folders
#

string(REGEX REPLACE "/+$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# folder for binaries and wrapper scripts
set(ESPRESSO_INSTALL_BINDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")

# folder for C++ and CUDA shared objects
set(ESPRESSO_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# python site-packages, can be overriden with CMake options
if(ESPRESSO_BUILD_WITH_PYTHON)
  if(NOT ESPRESSO_INSTALL_PYTHON)
    if(CMAKE_INSTALL_PREFIX STREQUAL "/")
      set(ESPRESSO_INSTALL_PYTHON "${Python_SITEARCH}")
    else()
      set(ESPRESSO_INSTALL_PYTHON
          "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages"
      )
    endif()
  endif()
  # override: package C++, CUDA and Cython shared objects together
  set(ESPRESSO_INSTALL_LIBDIR "${ESPRESSO_INSTALL_PYTHON}/espressomd")
endif()

#
# Libraries
#

if(ESPRESSO_BUILD_WITH_FFTW)
  find_package(FFTW3 REQUIRED)
endif()

# We need the parallel hdf5 version!
if(ESPRESSO_BUILD_WITH_HDF5)
  # The FindHDF5 function will fall back to the serial version if no parallel
  # version was found, and print to the CMake log that HDF5 was found. There is
  # no QUIET argument to override that message. This can be confusing to people
  # who are not familiar with the way hdf5 is distributed in Linux package
  # repositories (libhdf5-dev is the serial version).
  set(HDF5_PREFER_PARALLEL 1)
  find_package(HDF5 "1.8" REQUIRED COMPONENTS C)
  if(HDF5_FOUND)
    if(HDF5_IS_PARALLEL)
      add_feature_info(HDF5 ON "parallel")
    else()
      set(HDF5_FOUND FALSE)
      message(FATAL_ERROR "HDF5 parallel version not found.")
    endif()
  endif()
endif()

# Check for the h5xx submodule and try to check it out if not found or update it
# if found.
if(ESPRESSO_BUILD_WITH_HDF5 AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  # Try to find git
  find_package(Git)
  if(GIT_FOUND)
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/libs/h5xx/.git")
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --
                              libs/h5xx WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    else()
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update -- libs/h5xx
                      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
  endif()
endif()

if(ESPRESSO_BUILD_WITH_SCAFACOS)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SCAFACOS scafacos REQUIRED)
endif()

if(ESPRESSO_BUILD_WITH_GSL)
  find_package(GSL REQUIRED)
endif()

if(ESPRESSO_BUILD_WITH_STOKESIAN_DYNAMICS)
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
  FetchContent_GetProperties(stokesian_dynamics)
  set(STOKESIAN_DYNAMICS 1)
  if(NOT stokesian_dynamics_POPULATED)
    FetchContent_MakeAvailable(stokesian_dynamics)
  endif()
endif()

if(ESPRESSO_BUILD_WITH_VALGRIND)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(VALGRIND valgrind REQUIRED)
  if(VALGRIND_FOUND)
    message(STATUS "Found valgrind: ${VALGRIND_INCLUDE_DIRS}")
  endif()
endif()

#
# MPI
#

find_package(MPI 3.0 REQUIRED)
include(espresso_get_mpiexec_vendor)
espresso_get_mpiexec_vendor()
if(${ESPRESSO_MPIEXEC_VERSION} VERSION_GREATER_EQUAL
   ${ESPRESSO_MINIMAL_MPIEXEC_VERSION})
  message(
    STATUS
      "Found ${ESPRESSO_MPIEXEC_VENDOR}: ${MPIEXEC} (found suitable version \"${ESPRESSO_MPIEXEC_VERSION}\", minimum required is \"${ESPRESSO_MINIMAL_MPIEXEC_VERSION}\")"
  )
else()
  message(
    FATAL_ERROR
      "Could not find a suitable ${ESPRESSO_MPIEXEC_VENDOR} implementation (found unsuitable version \"${ESPRESSO_MPIEXEC_VERSION}\", minimum required is \"${ESPRESSO_MINIMAL_MPIEXEC_VERSION}\")"
  )
endif()

# OpenMPI checks the number of processes against the number of CPUs
set(ESPRESSO_MPIEXEC_OVERSUBSCRIBE "")
# Open MPI 4.x has a bug on NUMA archs that prevents running in singleton mode
set(ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA OFF)
set(ESPRESSO_CPU_MODEL_NAME_OMPI_SINGLETON_NUMA_PATTERN "AMD (EPYC|Ryzen)")

if("${ESPRESSO_MPIEXEC_VENDOR}" STREQUAL "OpenMPI")
  set(ESPRESSO_MPIEXEC_OVERSUBSCRIBE "-oversubscribe")
  if(${ESPRESSO_MPIEXEC_VERSION} VERSION_LESS 5.0)
    if(NOT DEFINED ESPRESSO_CPU_MODEL_NAME)
      if(CMAKE_SYSTEM_NAME STREQUAL Linux)
        if(EXISTS /proc/cpuinfo)
          file(READ /proc/cpuinfo ESPRESSO_CPU_INFO)
          string(REGEX
                 REPLACE ".*\n[Mm]odel name[ \t]*:[ \t]+([^\n]+).*" "\\1"
                         ESPRESSO_CPU_MODEL_NAME_STRING "${ESPRESSO_CPU_INFO}")
        else()
          set(ESPRESSO_CPU_MODEL_NAME_STRING "__unreadable")
        endif()
      else()
        set(ESPRESSO_CPU_MODEL_NAME_STRING "__unaffected")
      endif()
      set(ESPRESSO_CPU_MODEL_NAME "${ESPRESSO_CPU_MODEL_NAME_STRING}"
          CACHE INTERNAL "")
    endif()
    if(ESPRESSO_CPU_MODEL_NAME MATCHES
       "^${ESPRESSO_CPU_MODEL_NAME_OMPI_SINGLETON_NUMA_PATTERN}")
      set(ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA ON)
    endif()
  endif()
endif()

# OpenMPI cannot run two jobs in parallel in a Docker container, because the
# same base folder is used to store the process ids of multiple jobs. Since the
# base folder is deleted upon completion of a job, other jobs will fail when
# attempting to create subdirectories in the base folder.
# https://github.com/open-mpi/ompi/issues/8510
if("${ESPRESSO_MPIEXEC_VENDOR}" STREQUAL "OpenMPI" AND ESPRESSO_INSIDE_DOCKER)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
  function(espresso_set_mpiexec_tmpdir)
    set(ESPRESSO_MPIEXEC_TMPDIR --mca orte_tmpdir_base
                                "/tmp/ompi.${hostname}.$ENV{USER}.${ARGV0}"
        PARENT_SCOPE)
  endfunction()
else()
  function(espresso_set_mpiexec_tmpdir)
    set(ESPRESSO_MPIEXEC_TMPDIR "" PARENT_SCOPE)
  endfunction()
endif()

# mpiexec argument to disable pinning of processes to cores. Used in test suite
# as several tests are run at the same time.
if(NOT APPLE)
  set(ESPRESSO_MPIEXEC_DONT_PIN --bind-to none)
else()
  # not supported on Apple
  set(ESPRESSO_MPIEXEC_DONT_PIN "")
endif()

#
# Boost
#

list(APPEND BOOST_COMPONENTS mpi serialization filesystem system)

if(ESPRESSO_BUILD_TESTS)
  list(APPEND BOOST_COMPONENTS unit_test_framework)
endif()

find_package(Boost 1.74.0 REQUIRED ${BOOST_COMPONENTS})

#
# Paths
#

set(CMAKE_INSTALL_RPATH "${ESPRESSO_INSTALL_LIBDIR}")

#
# Flags
#

# drop 'lib' prefix from all libraries
set(CMAKE_SHARED_LIBRARY_PREFIX "")

add_library(espresso_coverage_flags INTERFACE)
add_library(espresso::coverage_flags ALIAS espresso_coverage_flags)
if(ESPRESSO_BUILD_WITH_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
      espresso_coverage_flags INTERFACE -g -fprofile-instr-generate
                                        -fcoverage-mapping)
  else()
    target_compile_options(espresso_coverage_flags INTERFACE -g --coverage
                                                             -fprofile-abs-path)
    target_link_libraries(espresso_coverage_flags INTERFACE gcov)
  endif()
endif()

add_library(espresso_cpp_flags INTERFACE)
add_library(espresso::cpp_flags ALIAS espresso_cpp_flags)
target_compile_options(
  espresso_cpp_flags
  INTERFACE
    -Wall
    -Wextra
    -pedantic
    # add extra warnings
    $<$<CXX_COMPILER_ID:Clang>:-Wextern-initializer>
    $<$<CXX_COMPILER_ID:Clang>:-Wrange-loop-analysis>
    -Wfloat-conversion
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wimplicit-float-conversion>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wunused-exception-parameter>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wmissing-variable-declarations>
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11.0.0>>:-Wnon-c-typedef-for-linkage>
    $<$<NOT:$<CXX_COMPILER_ID:Intel>>:-Wdelete-non-virtual-dtor>
    # disable warnings from -Wall and -Wextra
    -Wno-sign-compare
    -Wno-unused-function
    -Wno-unused-parameter
    -Wno-array-bounds
    $<$<CXX_COMPILER_ID:GNU>:-Wno-restrict>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-clobbered>
    $<$<CXX_COMPILER_ID:Intel,IntelLLVM>:-diag-disable=592>
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-gnu-zero-variadic-macro-arguments>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,8.1.0>>:-Wno-cast-function-type>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-implicit-fallthrough>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM,GNU>>:-Wno-unused-private-field>
    # warnings are errors
    $<$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>:-Werror>)

if(ESPRESSO_BUILD_WITH_CUDA)
  target_compile_options(
    espresso_cuda_flags
    INTERFACE
      -Wall
      -Wextra
      -Wno-sign-compare
      -Wno-unused-parameter
      $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-implicit-fallthrough>
      # warnings are errors
      $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<CUDA_COMPILER_ID:NVIDIA>>:--Werror=all-warnings>
      $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<CUDA_COMPILER_ID:Clang>>:-Werror>
  )
endif()

# disable warning from -Wextra on ARM processors
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_PROCESSOR MATCHES
                                            "arm")
  target_compile_options(espresso_cpp_flags INTERFACE -Wno-psabi)
endif()

# prevent 80-bit arithmetic on old Intel processors
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SIZEOF_VOID_P EQUAL 4
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "[xX]86")
  target_compile_options(espresso_cpp_flags INTERFACE -ffloat-store)
endif()

# enable boost::variant with more than 20 types
target_compile_options(
  espresso_cpp_flags INTERFACE -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
                               -DBOOST_MPL_LIMIT_LIST_SIZE=30)

set(CMAKE_MACOSX_RPATH TRUE)

#
# Sanitizers
#

if(ESPRESSO_BUILD_WITH_ASAN AND ESPRESSO_BUILD_WITH_MSAN)
  message(
    FATAL_ERROR
      "Address sanitizer and memory sanitizer cannot be enabled simultaneously")
endif()
if(ESPRESSO_BUILD_WITH_ASAN)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O1")
  target_compile_options(espresso_cpp_flags INTERFACE -fsanitize=address
                                                      -fno-omit-frame-pointer)
  target_link_libraries(espresso_cpp_flags INTERFACE -fsanitize=address)
  if(ESPRESSO_BUILD_WITH_CUDA)
    target_link_libraries(espresso_cuda_flags INTERFACE -fsanitize=address)
  endif()
endif()
if(ESPRESSO_BUILD_WITH_MSAN)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O1")
  target_compile_options(espresso_cpp_flags INTERFACE -fsanitize=memory
                                                      -fno-omit-frame-pointer)
  target_link_libraries(espresso_cpp_flags INTERFACE -fsanitize=memory)
  if(ESPRESSO_BUILD_WITH_CUDA)
    target_link_libraries(espresso_cuda_flags INTERFACE -fsanitize=memory)
  endif()
endif()
if(ESPRESSO_BUILD_WITH_UBSAN)
  target_compile_options(espresso_cpp_flags INTERFACE -fsanitize=undefined)
  target_link_libraries(espresso_cpp_flags INTERFACE -fsanitize=undefined)
  if(ESPRESSO_BUILD_WITH_CUDA)
    target_link_libraries(espresso_cuda_flags INTERFACE -fsanitize=undefined)
  endif()
endif()

target_link_libraries(espresso_cpp_flags INTERFACE espresso::coverage_flags)

#
# Static analysis
#

if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
  find_package(ClangTidy "${CMAKE_CXX_COMPILER_VERSION}" EXACT REQUIRED)
  set(ESPRESSO_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  set(ESPRESSO_CUDA_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=--cuda-host-only")
endif()

if(ESPRESSO_BUILD_WITH_CPPCHECK)
  find_program(CMAKE_CXX_CPPCHECK NAMES cppcheck)
  if(NOT CMAKE_CXX_CPPCHECK)
    message(FATAL_ERROR "Could not find the program cppcheck.")
  endif()
  list(APPEND CMAKE_CXX_CPPCHECK "--enable=all"
       "--std=c++${CMAKE_CXX_STANDARD}" "--quiet" "--inline-suppr"
       "--suppressions-list=${CMAKE_CURRENT_SOURCE_DIR}/.cppcheck")
  if(ESPRESSO_WARNINGS_ARE_ERRORS)
    list(APPEND CMAKE_CXX_CPPCHECK "--error-exitcode=2")
  endif()
endif()

#
# Testing
#

if(ESPRESSO_BUILD_TESTS)
  enable_testing()
  add_custom_target(check)
  set(ESPRESSO_CTEST_ARGS ""
      CACHE STRING
            "Extra arguments to give to ctest calls (separated by semicolons)")
  set(ESPRESSO_TEST_NP "4" CACHE STRING
                                 "Maximal number of MPI ranks to use per test")
  if(ESPRESSO_BUILD_WITH_PYTHON)
    add_subdirectory(testsuite)
  endif()
endif()

if(ESPRESSO_BUILD_BENCHMARKS)
  add_custom_target(benchmark)
  add_subdirectory(maintainer/benchmarks)
endif()

#
# waLBerla
#

if(ESPRESSO_BUILD_WITH_WALBERLA)
  FetchContent_GetProperties(walberla)
  string(REGEX REPLACE "([/\\]walberla)-src$" "\\1-build" walberla_BINARY_DIR
                       "${walberla_SOURCE_DIR}")
  set(WALBERLA_BUILD_TESTS off CACHE BOOL "")
  set(WALBERLA_BUILD_BENCHMARKS off CACHE BOOL "")
  set(WALBERLA_BUILD_TOOLS off CACHE BOOL "")
  set(WALBERLA_BUILD_TUTORIALS off CACHE BOOL "")
  set(WALBERLA_BUILD_SHOWCASES off CACHE BOOL "")
  set(WALBERLA_BUILD_DOC off CACHE BOOL "")
  set(WALBERLA_LOGLEVEL "WARNING" CACHE STRING "")
  set(CMAKE_POSITION_INDEPENDENT_CODE on CACHE BOOL "")
  if(ESPRESSO_BUILD_WITH_CUDA)
    set(WALBERLA_BUILD_WITH_CUDA "on" CACHE BOOL "")
    if(NOT ESPRESSO_CUDA_COMPILER STREQUAL "clang")
      if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "variable CMAKE_CUDA_ARCHITECTURES is undefined")
      endif()
    endif()
  endif()
  if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
    set(ESPRESSO_USE_WALBERLA_FFT 1)
    set(WALBERLA_BUILD_WITH_FFTW on CACHE BOOL "")
  else()
    set(WALBERLA_BUILD_WITH_FFTW off CACHE BOOL "")
  endif()
  set(WALBERLA_BUILD_WITH_FASTMATH off CACHE BOOL "")
  if(NOT walberla_POPULATED)
    FetchContent_MakeAvailable(walberla)
  endif()
  set(WALBERLA_LIBS
      walberla::core walberla::domain_decomposition walberla::blockforest
      walberla::boundary walberla::field walberla::lbm walberla::timeloop
      walberla::vtk)
  if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
    set(WALBERLA_LIBS ${WALBERLA_LIBS} walberla::fft)
  endif()
  if(ESPRESSO_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
    set(WALBERLA_LIBS ${WALBERLA_LIBS} walberla::gpu)
  endif()
  # workaround for https://gitlab.kitware.com/cmake/cmake/-/issues/21283
  foreach(target_w_namespace IN LISTS WALBERLA_LIBS)
    string(REPLACE "walberla::" "" target_wo_namespace ${target_w_namespace})
    add_library(${target_w_namespace} ALIAS ${target_wo_namespace})
  endforeach()
  if(ESPRESSO_BUILD_WITH_WALBERLA_AVX)
    function(espresso_avx_flags_callback COMPILER_AVX2_FLAG)
      target_compile_options(
        espresso_avx_flags INTERFACE "${COMPILER_AVX2_FLAG}"
                                     "-DESPRESSO_BUILD_WITH_AVX_KERNELS")
    endfunction()
    espresso_enable_avx2_support(espresso_avx_flags_callback)
  endif()
endif()

if(ESPRESSO_BUILD_WITH_CALIPER)
  FetchContent_GetProperties(caliper)
  set(CALIPER_OPTION_PREFIX on CACHE BOOL "")
  set(CALIPER_WITH_MPI on CACHE BOOL "")
  set(CALIPER_WITH_NVTX off CACHE BOOL "")
  set(CALIPER_WITH_CUPTI off CACHE BOOL "")
  set(CALIPER_BUILD_SHARED_LIBS on CACHE BOOL "")
  if(NOT caliper_POPULATED)
    FetchContent_MakeAvailable(caliper)
  endif()
  target_compile_options(
    caliper-services
    PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  target_compile_options(
    caliper-runtime
    PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
            $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-volatile>)
endif()

#
# Subdirectories
#

add_subdirectory(doc)
add_subdirectory(src)
add_subdirectory(libs)

#
# Feature summary
#

include(FeatureSummary)
feature_summary(WHAT ALL)
