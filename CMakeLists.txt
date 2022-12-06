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

cmake_minimum_required(VERSION 3.17)
message(STATUS "CMake version: ${CMAKE_VERSION}")
if(POLICY CMP0076)
  # make target_sources() convert relative paths to absolute
  cmake_policy(SET CMP0076 NEW)
endif()
if(POLICY CMP0025)
  # make CXX_COMPILER_ID return "AppleClang" instead of "Clang" for Apple Clang
  cmake_policy(SET CMP0025 NEW)
endif()
# CMake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# C++ standard
enable_language(CXX)
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# require minimal compiler versions
macro(minimal_compiler_version)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "${ARGV0}" AND CMAKE_CXX_COMPILER_VERSION
                                                   VERSION_LESS "${ARGV1}")
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID} \
                 ${CMAKE_CXX_COMPILER_VERSION} (required version >= ${ARGV1})")
  endif()
endmacro()

minimal_compiler_version("GNU" 7.0.0)
minimal_compiler_version("Clang" 8.0.0)
minimal_compiler_version("AppleClang" 11.0.0)
minimal_compiler_version("Intel" 18.0)
minimal_compiler_version("IntelLLVM" 2021.0)

include(FeatureSummary)
project(ESPResSo)
include(GNUInstallDirs)
include(option_enum)
if(POLICY CMP0074)
  # make find_package() use <PackageName>_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif()

set(PROJECT_VERSION "4.3-dev")

#
# CMake internal vars
#

# Select the build type
option_enum(
  varname "CMAKE_BUILD_TYPE" help_text "build type" default_value "Release"
  possible_values
  "Debug;Release;RelWithDebInfo;MinSizeRel;Coverage;RelWithAssert")
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE} -Og -g")
set(CMAKE_CXX_FLAGS_RELWITHASSERT "${CMAKE_CXX_FLAGS_RELWITHASSERT} -O3 -g")

# On Mac OS X, first look for other packages, then frameworks
set(CMAKE_FIND_FRAMEWORK LAST)

# ##############################################################################
# User input options
# ##############################################################################

option(WITH_PYTHON "Build with Python bindings" ON)
option(WITH_GSL "Build with GSL support" OFF)
option(WITH_CUDA "Build with GPU support" OFF)
option(WITH_HDF5 "Build with HDF5 support" OFF)
option(WITH_TESTS "Enable tests" ON)
option(WITH_SCAFACOS "Build with ScaFaCoS support" OFF)
option(WITH_STOKESIAN_DYNAMICS "Build with Stokesian Dynamics" OFF)
option(WITH_WALBERLA "Build with WALBERLA Lattice-Boltzmann support" OFF)
option(WALBERLA_USE_AVX
       "Build WALBERLA Lattice-Boltzmann with AVX vectorization" OFF)
option(ESPRESSO_BUILD_WITH_WALBERLA_FFT "Build waLBerla with FFT support" ${WITH_WALBERLA})
option(WITH_BENCHMARKS "Enable benchmarks" OFF)
option(WITH_VALGRIND_INSTRUMENTATION
       "Build with valgrind instrumentation markers" OFF)
option(WITH_CPPCHECK "Run Cppcheck during compilation" OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  option(WITH_CLANG_TIDY "Run Clang-Tidy during compilation" OFF)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                            "GNU")
  option(WITH_COVERAGE "Generate code coverage report for C++ code" OFF)
  option(WITH_COVERAGE_PYTHON "Generate code coverage report for Python code"
         OFF)
  option(WITH_ASAN "Build with address sanitizer" OFF)
  option(WITH_UBSAN "Build with undefined behavior sanitizer" OFF)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT APPLE)
  option(
    WITH_MSAN
    "Build with memory sanitizer (experimental; requires a memory-sanitized Python interpreter)"
    OFF)
endif()
option(
  ESPRESSO_ADD_OMPI_SINGLETON_WARNING
  "Add a runtime warning in the pypresso script for NUMA architectures that aren't supported in singleton mode by Open MPI 4.x"
  ON)
option(WARNINGS_ARE_ERRORS "Treat warnings as errors during compilation" OFF)
option(WITH_CCACHE "Use ccache compiler invocation." OFF)
option(WITH_PROFILER "Enable profiler annotations." OFF)
option(INSIDE_DOCKER "Enable when running inside Docker." OFF)
set(TEST_TIMEOUT "300" CACHE STRING
                             "Timeout in seconds for each testsuite test")

if(WITH_CCACHE)
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
set(MYCONFIG_NAME "myconfig.hpp" CACHE STRING
                                       "Default name of the local config file")

# Check which config file to use
include(MyConfig)

if(WITH_WALBERLA)
  set(CMAKE_CXX_STANDARD 17)
endif()

#
# Pretty function
#

include(CheckCXXSourceCompiles)

# cross-platform macro to print the function name in error messages
set(PRETTY_FUNCTION_EXTENSION __func__)

# search for a supported compiler extension that prints the function name as
# well as its list of arguments, return type and namespace
foreach(func_name __PRETTY_FUNCTION__ __FUNCSIG__ __FUNCTION__)
  check_cxx_source_compiles(
    "
     #include <string>
     int main() { std::string(${func_name}); }
     " result${func_name})
  if(result${func_name})
    set(PRETTY_FUNCTION_EXTENSION ${func_name})
    break()
  endif(result${func_name})
endforeach()

#
# AVX2 support
#

include(CheckCXXCompilerFlag)

add_library(espresso_avx_flags INTERFACE)
add_library(espresso::avx_flags ALIAS espresso_avx_flags)

function(enable_avx2_support)
  set(compiler_avx2_flag "")
  foreach(flag_name "-mavx2" "/arch:AVX2")
    string(REGEX REPLACE "[^0-9A-Za-z_]" "_" flag_variable "${flag_name}")
    check_cxx_compiler_flag("${flag_name}"
                            compiler_has_${flag_variable}_flag_result)
    if(compiler_has_${flag_variable}_flag_result)
      set(compiler_avx2_flag "${flag_name}")
      target_compile_options(espresso_avx_flags
                             INTERFACE "${compiler_avx2_flag}" "-DWITH_AVX_KERNELS")
      break()
    endif()
  endforeach()
  if(compiler_avx2_flag STREQUAL "")
    message(
      FATAL_ERROR
        "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} doesn't support AVX2-specific compiler flags."
    )
  endif()
  if(NOT compiler_avx2_flag STREQUAL "/arch:AVX2")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -march=native -E -v - INPUT_FILE /dev/null
      OUTPUT_VARIABLE march_native_output_string
      ERROR_VARIABLE march_native_output_string)
    if(NOT "${march_native_output_string}" MATCHES
       "[ \n](\\+avx2|-mavx2|-D__AVX2__)[ \n]")
      message(
        FATAL_ERROR
          "AVX2 not supported on this CPU architecture according to ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. While ESPResSo will still compile, you will get SIGILL crashes when calling AVX functions."
      )
    endif()
  endif()
  set(CMAKE_REQUIRED_FLAGS_BACKUP "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${compiler_avx2_flag}")
  check_cxx_source_compiles(
    "#include <immintrin.h>
       __m256i xi_i = _mm256_set_epi32(1, 2, 3, 4, 5, 6, 7, 8);
       __m256  xi_s = _mm256_set_ps(0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f);
       __m256d xi_d = _mm256_set_pd(0.0, 1.0, 2.0, 3.0);
       int main() {}
      " compiler_has_avx2_support)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_BACKUP}")
  if(NOT compiler_has_avx2_support)
    message(
      FATAL_ERROR
        "Cannot execute a simple AVX2 program compiled by ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}."
    )
  endif()
endfunction()

#
# Interface libraries
#

# CUDA compiler
if(WITH_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(MINIMAL_CUDA_VERSION 10.0)
  find_package(CUDAToolkit ${MINIMAL_CUDA_VERSION} REQUIRED)
  option_enum(varname "WITH_CUDA_COMPILER" help_text "CUDA compiler"
              default_value "nvcc" possible_values "nvcc;clang")
  if(WITH_CUDA_COMPILER STREQUAL "nvcc")
    find_package(CUDACompilerNVCC ${MINIMAL_CUDA_VERSION} REQUIRED)
  elseif(WITH_CUDA_COMPILER STREQUAL "clang")
    if(WITH_COVERAGE)
      message(
        FATAL_ERROR
          "Cannot enable code coverage with Clang as the CUDA compiler")
    endif()
    find_package(CUDACompilerClang 9.0 REQUIRED)
  else()
    message(FATAL_ERROR "Unknown CUDA compiler '${WITH_CUDA_COMPILER}'")
  endif()
endif(WITH_CUDA)

# Python interpreter and Cython interface library
if(WITH_PYTHON)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Cython 0.29.14 REQUIRED)
  find_program(IPYTHON_EXECUTABLE NAMES jupyter ipython3 ipython)
endif(WITH_PYTHON)

#
# Installation folders
#

string(REGEX REPLACE "/+$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# folder for binaries and wrapper scripts
set(ESPRESSO_INSTALL_BINDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")

# folder for C++ and CUDA shared objects
set(ESPRESSO_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# python site-packages, can be overriden with CMake options
if(WITH_PYTHON)
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

find_package(FFTW3)
if(FFTW3_FOUND)
  set(FFTW 3)
endif(FFTW3_FOUND)

# We need the parallel hdf5 version!
if(WITH_HDF5)
  # The FindHDF5 function will fall back to the serial version if no parallel
  # version was found, and print to the CMake log that HDF5 was found. There is
  # no QUIET argument to override that message. This can be confusing to people
  # who are not familiar with the way hdf5 is distributed in Linux package
  # repositories (libhdf5-dev is the serial version).
  set(HDF5_PREFER_PARALLEL 1)
  find_package(HDF5 "1.8" REQUIRED COMPONENTS C)
  if(HDF5_FOUND)
    if(HDF5_IS_PARALLEL)
      set(H5MD 1)
      add_feature_info(HDF5 ON "parallel")
    else()
      unset(H5MD)
      set(HDF5_FOUND FALSE)
      message(FATAL_ERROR "HDF5 parallel version not found.")
    endif(HDF5_IS_PARALLEL)
  endif(HDF5_FOUND)
endif(WITH_HDF5)

# Check for the h5xx submodule and try to check it out if not found or update it
# if found.
if(WITH_HDF5 AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
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

if(WITH_SCAFACOS)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SCAFACOS scafacos REQUIRED)
  if(SCAFACOS_FOUND)
    set(SCAFACOS 1)
  endif(SCAFACOS_FOUND)
endif(WITH_SCAFACOS)

if(WITH_GSL)
  find_package(GSL REQUIRED)
else()
  find_package(GSL)
endif(WITH_GSL)

if(GSL_FOUND)
  set(GSL 1)
endif(GSL_FOUND)

if(WITH_STOKESIAN_DYNAMICS)
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
  include(FetchContent)
  FetchContent_Declare(
    stokesian_dynamics
    GIT_REPOSITORY https://github.com/hmenke/espresso-stokesian-dynamics.git
    GIT_TAG 862a7537a366f0c32f0c25e46bd107bea590faea)
  FetchContent_GetProperties(stokesian_dynamics)
  set(STOKESIAN_DYNAMICS 1)
  if(NOT stokesian_dynamics_POPULATED)
    FetchContent_Populate(stokesian_dynamics)
    add_subdirectory(${stokesian_dynamics_SOURCE_DIR}
                     ${stokesian_dynamics_BINARY_DIR})
  endif()
endif(WITH_STOKESIAN_DYNAMICS)

if(WITH_VALGRIND_INSTRUMENTATION)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(VALGRIND valgrind REQUIRED)
  if(VALGRIND_FOUND)
    set(VALGRIND_INSTRUMENTATION 1)
    message(STATUS ${VALGRIND_INCLUDE_DIRS})
    include_directories(SYSTEM ${VALGRIND_INCLUDE_DIRS})
  endif(VALGRIND_FOUND)
endif(WITH_VALGRIND_INSTRUMENTATION)

#
# MPI
#

find_package(MPI 3.0 REQUIRED)
include(espresso_get_mpiexec_vendor)
espresso_get_mpiexec_vendor()
message(
  STATUS
    "Found ${ESPRESSO_MPIEXEC_VENDOR}: ${MPIEXEC} (version \"${ESPRESSO_MPIEXEC_VERSION}\")"
)

# OpenMPI checks the number of processes against the number of CPUs
set(MPIEXEC_OVERSUBSCRIBE "")
# Open MPI 4.x has a bug on NUMA archs that prevents running in singleton mode
set(ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA OFF)
set(ESPRESSO_CPU_MODEL_NAME_OMPI_SINGLETON_NUMA_PATTERN "AMD (EPYC|Ryzen)")

if("${ESPRESSO_MPIEXEC_VENDOR}" STREQUAL "OpenMPI")
  if(${ESPRESSO_MPIEXEC_VERSION} VERSION_GREATER_EQUAL 2.0.0)
    set(MPIEXEC_OVERSUBSCRIBE "-oversubscribe")
  endif()
  if(${ESPRESSO_MPIEXEC_VERSION} VERSION_GREATER_EQUAL 4.0
     AND ${ESPRESSO_MPIEXEC_VERSION} VERSION_LESS 5.0)
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
if("${ESPRESSO_MPIEXEC_VENDOR}" STREQUAL "OpenMPI" AND INSIDE_DOCKER)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
  function(set_mpiexec_tmpdir)
    set(MPIEXEC_TMPDIR --mca orte_tmpdir_base
                       "/tmp/ompi.${hostname}.$ENV{USER}.${ARGV0}" PARENT_SCOPE)
  endfunction()
else()
  function(set_mpiexec_tmpdir)
    set(MPIEXEC_TMPDIR "" PARENT_SCOPE)
  endfunction()
endif()

#
# Boost
#

list(APPEND BOOST_COMPONENTS mpi serialization filesystem system)

if(WITH_TESTS)
  list(APPEND BOOST_COMPONENTS unit_test_framework)
endif()

find_package(Boost 1.69.0 REQUIRED ${BOOST_COMPONENTS})

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
if(WITH_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
      espresso_coverage_flags INTERFACE -g -fprofile-instr-generate
                                        -fcoverage-mapping)
  else()
    target_compile_options(
      espresso_coverage_flags INTERFACE -g --coverage -fprofile-arcs
                                        -ftest-coverage)
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
    $<$<BOOL:${WARNINGS_ARE_ERRORS}>:-Werror>
    # add extra warnings
    $<$<CXX_COMPILER_ID:Clang>:-Wextern-initializer>
    $<$<CXX_COMPILER_ID:Clang>:-Wrange-loop-analysis>
    -Wfloat-conversion
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wimplicit-float-conversion>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wunused-exception-parameter>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wmissing-variable-declarations>
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11.0.0>>:-Wnon-c-typedef-for-linkage>
    $<$<NOT:$<CXX_COMPILER_ID:Intel>>:-Wdelete-non-virtual-dtor>
    # disable warnings from -Wextra
    -Wno-sign-compare
    -Wno-unused-function
    -Wno-unused-parameter
    $<$<CXX_COMPILER_ID:GNU>:-Wno-clobbered>
    $<$<CXX_COMPILER_ID:Intel,IntelLLVM>:-diag-disable=592>
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-gnu-zero-variadic-macro-arguments>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,8.1.0>>:-Wno-cast-function-type>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-implicit-fallthrough>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM,GNU>>:-Wno-unused-private-field>)

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

if(WITH_ASAN AND WITH_MSAN)
  message(
    FATAL_ERROR
      "Address sanitizer and memory sanitizer cannot be enabled simultaneously")
endif()
if(WITH_ASAN)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O1")
  target_compile_options(espresso_cpp_flags INTERFACE -fsanitize=address
                                                      -fno-omit-frame-pointer)
  target_link_libraries(espresso_cpp_flags INTERFACE -fsanitize=address)
endif()
if(WITH_MSAN)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O1")
  target_compile_options(espresso_cpp_flags INTERFACE -fsanitize=memory
                                                      -fno-omit-frame-pointer)
  target_link_libraries(espresso_cpp_flags INTERFACE -fsanitize=memory)
endif()
if(WITH_UBSAN)
  target_compile_options(espresso_cpp_flags INTERFACE -fsanitize=undefined)
  target_link_libraries(espresso_cpp_flags INTERFACE -fsanitize=undefined)
endif()

target_link_libraries(espresso_cpp_flags INTERFACE espresso::coverage_flags)

#
# Static analysis
#

if(WITH_CLANG_TIDY)
  find_package(ClangTidy "${CMAKE_CXX_COMPILER_VERSION}" EXACT REQUIRED)
  set(ESPRESSO_DEFAULT_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=--cuda-host-only")
  set(CMAKE_CXX_CLANG_TIDY "${ESPRESSO_DEFAULT_CXX_CLANG_TIDY}")
endif()

if(WITH_CPPCHECK)
  find_program(CMAKE_CXX_CPPCHECK NAMES cppcheck)
  if(NOT CMAKE_CXX_CPPCHECK)
    message(FATAL_ERROR "Could not find the program cppcheck.")
  endif()
  list(APPEND CMAKE_CXX_CPPCHECK "--enable=all"
       "--std=c++${CMAKE_CXX_STANDARD}" "--quiet" "--inline-suppr"
       "--suppressions-list=${CMAKE_CURRENT_SOURCE_DIR}/.cppcheck")
  if(WARNINGS_ARE_ERRORS)
    list(APPEND CMAKE_CXX_CPPCHECK "--error-exitcode=2")
  endif()
endif()

#
# Testing
#

if(WITH_TESTS)
  enable_testing()
  add_custom_target(check)
  set(CTEST_ARGS ""
      CACHE STRING
            "Extra arguments to give to ctest calls (separated by semicolons)")
  set(TEST_NP "4" CACHE STRING "Maximal number of MPI ranks to use per test")
  if(WITH_PYTHON)
    add_subdirectory(testsuite)
  endif(WITH_PYTHON)
endif(WITH_TESTS)

if(WITH_BENCHMARKS)
  add_custom_target(benchmark)
  add_subdirectory(maintainer/benchmarks)
endif(WITH_BENCHMARKS)

# Walberla
if(WITH_WALBERLA)
  include(FetchContent)
  FetchContent_Declare(
    walberla
    GIT_REPOSITORY https://i10git.cs.fau.de/walberla/walberla.git
    GIT_TAG        01a28162ae1aacf7b96152c9f886ce54cc7f53ff
    PATCH_COMMAND  patch -p1 -sN --ignore-whitespace -i
                   ${PROJECT_SOURCE_DIR}/cmake/waLBerla.patch || echo "Patch not applied!"
  )
  if(NOT walberla_POPULATED)
    FetchContent_Populate(walberla)
    set(WALBERLA_BUILD_TESTS off CACHE BOOL "")
    set(WALBERLA_BUILD_BENCHMARKS off CACHE BOOL "")
    set(WALBERLA_BUILD_TOOLS off CACHE BOOL "")
    set(WALBERLA_BUILD_TUTORIALS off CACHE BOOL "")
    set(WALBERLA_BUILD_SHOWCASES off CACHE BOOL "")
    set(WALBERLA_BUILD_DOC off CACHE BOOL "")
    set(CMAKE_POSITION_INDEPENDENT_CODE on CACHE BOOL "")
    if (WITH_CUDA AND CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11)
      set(WALBERLA_BUILD_WITH_CUDA "on"  CACHE BOOL "")
#      if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
#        set(CMAKE_CUDA_ARCHITECTURES 75)
#      endif()
    endif()
    if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
      set(ESPRESSO_USE_WALBERLA_FFT 1)
      set(WALBERLA_BUILD_WITH_FFTW on CACHE BOOL "")
    else()
      set(WALBERLA_BUILD_WITH_FFTW off CACHE BOOL "")
    endif()
    set(CMAKE_CXX_CLANG_TIDY "")
    add_subdirectory("${walberla_SOURCE_DIR}" "${walberla_BINARY_DIR}")
    set(CMAKE_CXX_CLANG_TIDY "${ESPRESSO_DEFAULT_CXX_CLANG_TIDY}")
    set(LB_WALBERLA 1)
    set(WALBERLA_LIBS
        walberla::core walberla::domain_decomposition walberla::blockforest
        walberla::boundary walberla::core walberla::domain_decomposition
        walberla::field walberla::lbm walberla::timeloop walberla::vtk)
    if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
      set(WALBERLA_LIBS ${WALBERLA_LIBS} walberla::fft)
    endif()
    # workaround for https://gitlab.kitware.com/cmake/cmake/-/issues/21283
    foreach(target_w_namespace IN LISTS WALBERLA_LIBS)
      string(REPLACE "walberla::" "" target_wo_namespace ${target_w_namespace})
      add_library(${target_w_namespace} ALIAS ${target_wo_namespace})
    endforeach()
    # TODO WALBERLA: remove -lstdc++fs
    # TODO WALBERLA: remove -lmpi_cxx -lmpi (required to find symbol MPI::Datatype::Free() in VTK code)
    # TODO WALBERLA: remove -lpfft (required to find pfft_free in FFT-code)
    # TODO WALBERLA: remove -lfftw3_mpi (required to find fftw_mpi_execute_r2r in FFT-code)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
      set(WALBERLA_LIBS ${WALBERLA_LIBS} stdc++fs)
      if(WALBERLA_BUILD_WITH_MPI)
        set(WALBERLA_LIBS ${WALBERLA_LIBS} -lmpi)
        if("${MPIEXEC_BACKEND_NAME}" STREQUAL "OpenMPI")
          set(WALBERLA_LIBS ${WALBERLA_LIBS} -lmpi_cxx)
        endif()
        if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
          set(WALBERLA_LIBS ${WALBERLA_LIBS} -lpfft -lfftw3_mpi)
        endif()
      endif()
    endif()
  endif()
  if(WALBERLA_USE_AVX)
    enable_avx2_support()
  endif()
endif(WITH_WALBERLA)

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
