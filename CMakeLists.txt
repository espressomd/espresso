# Copyright (C) 2009-2019 The ESPResSo project
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

cmake_minimum_required(VERSION 3.16)
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

include(FeatureSummary)
include(GNUInstallDirs)
project(ESPResSo)
include(option_enum)
if(POLICY CMP0074)
  # make find_package() use <PackageName>_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif()

set(PROJECT_VERSION "4.2-dev")

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
option(WITH_BENCHMARKS "Enable benchmarks" OFF)
option(WITH_VALGRIND_INSTRUMENTATION
       "Build with valgrind instrumentation markers" OFF)
option(WITH_WALBERLA "Build with WALBERLA Lattice-Boltzmann support" OFF)
option(WALBERLA_USE_AVX "Build WALBERLA Lattice-Boltzmann with AVX vectorization" OFF)
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
option(WARNINGS_ARE_ERRORS "Treat warnings as errors during compilation" OFF)
option(WITH_CCACHE "Use ccache compiler invocation." OFF)
option(WITH_PROFILER "Enable profiler annotations." OFF)
set(TEST_TIMEOUT "300" CACHE STRING
                             "Timeout in seconds for each testsuite test")

if(WITH_CCACHE)
  find_program(CCACHE ccache)
  if(CCACHE)
    message(STATUS "Found ccache ${CCACHE}")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE} CACHE STRING "ccache executable")
  else()
    message(FATAL_ERROR "ccache not found.")
  endif(CCACHE)
endif(WITH_CCACHE)

# Write compile commands to file, for various tools...
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# choose the name of the config file
set(MYCONFIG_NAME "myconfig.hpp" CACHE STRING
                                       "Default name of the local config file")

# Check which config file to use
include(MyConfig)

#
# Pretty function
#

include(CheckCXXSourceCompiles)

set(__PRETTYFUNC__ __func__)
foreach(func_name __PRETTY_FUNCTION__ __FUNCTION__)
  check_cxx_source_compiles(
    "
     #include <string>
     int main() { std::string(${func_name}); }
     " result${func_name})
  if(result${func_name})
    set(__PRETTYFUNC__ ${func_name})
    break()
  endif(result${func_name})
endforeach(func_name __PRETTY_FUNCTION__ __FUNCTION__)

#
# Libraries
#

# CUDA compiler
if(WITH_CUDA)
  set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
  set(MINIMAL_CUDA_VERSION 10.0)
  option_enum(varname "WITH_CUDA_COMPILER" help_text "CUDA compiler"
              default_value "nvcc" possible_values "nvcc;clang")
  if(WITH_CUDA_COMPILER STREQUAL "nvcc")
    find_package(CUDA ${MINIMAL_CUDA_VERSION} REQUIRED)
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

find_package(PythonInterp 3.7 REQUIRED)

if(WITH_PYTHON)
  find_package(Cython 0.29.14 REQUIRED)
  find_package(PythonHeaders REQUIRED)
  find_package(NumPy REQUIRED)
  find_program(IPYTHON_EXECUTABLE NAMES jupyter ipython3 ipython)
endif(WITH_PYTHON)

find_package(FFTW3)
if(FFTW3_FOUND)
  set(FFTW 3)
endif(FFTW3_FOUND)

# If we build Python bindings, turn on script interface
if(WITH_PYTHON)
  set(WITH_SCRIPT_INTERFACE ON)
endif()

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
  set(CMAKE_INSTALL_LIBDIR
      "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR}/espressomd")
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

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR}/espressomd")

#
# Flags
#

# drop 'lib' prefix from all libraries
set(CMAKE_SHARED_LIBRARY_PREFIX "")

add_library(coverage_interface INTERFACE)
if(WITH_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
      coverage_interface INTERFACE -g -fprofile-instr-generate
                                   -fcoverage-mapping)
  else()
    target_compile_options(
      coverage_interface INTERFACE -g --coverage -fprofile-arcs -ftest-coverage)
    target_link_libraries(coverage_interface INTERFACE gcov)
  endif()
endif()

add_library(cxx_interface INTERFACE)
target_compile_options(
  cxx_interface
  INTERFACE -Wall
            -Wextra
            $<$<BOOL:${WARNINGS_ARE_ERRORS}>:-Werror>
            # add extra warnings
            $<$<CXX_COMPILER_ID:Clang>:-Wextern-initializer>
            $<$<CXX_COMPILER_ID:Clang>:-Wrange-loop-analysis>
            -Wfloat-conversion
            # disable warnings from -Wextra
            -Wno-sign-compare
            -Wno-unused-function
            -Wno-unused-variable
            -Wno-unused-parameter
            -Wno-missing-braces
            $<$<CXX_COMPILER_ID:GNU>:-Wno-clobbered>
            $<$<CXX_COMPILER_ID:Intel>:-wd592>)

# disable more warnings from -Wextra
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_PROCESSOR MATCHES
                                            "arm")
  target_compile_options(cxx_interface INTERFACE -Wno-psabi)
endif()
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel"
   AND NOT (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
                                                     VERSION_LESS "7.0.0"))
  target_compile_options(cxx_interface INTERFACE -Wno-implicit-fallthrough)
endif()
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_ID
                                                STREQUAL "Intel")
  target_compile_options(cxx_interface INTERFACE -Wno-unused-private-field)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(cxx_interface
                         INTERFACE -Wno-gnu-zero-variadic-macro-arguments)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "8.0.0")
    target_compile_options(cxx_interface INTERFACE -Wimplicit-float-conversion)
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
     AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "11.0.0")
    target_compile_options(cxx_interface INTERFACE -Wnon-c-typedef-for-linkage)
  endif()
endif()

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_VERSION
                                               VERSION_GREATER "4.8.5")
  # older GCC versions don't support -Wno-pedantic which we need in src/python
  target_compile_options(cxx_interface INTERFACE -pedantic)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
                                            VERSION_GREATER_EQUAL "8.1.0")
  target_compile_options(cxx_interface INTERFACE -Wno-cast-function-type)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND CMAKE_CXX_COMPILER_VERSION
                                              VERSION_LESS "16.0")
  # workaround for compiler crash related to decltype() and variadic template
  # usage inside Boost
  target_compile_options(cxx_interface
                         INTERFACE -DBOOST_NO_CXX11_VARIADIC_TEMPLATES)
endif()

# prevent 80-bit arithmetic on old Intel processors
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SIZEOF_VOID_P EQUAL 4
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "[xX]86")
  target_compile_options(cxx_interface INTERFACE -ffloat-store)
endif()

# enable boost::variant with more than 20 types
target_compile_options(
  cxx_interface INTERFACE -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
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
  target_compile_options(cxx_interface INTERFACE -fsanitize=address
                                                 -fno-omit-frame-pointer)
  target_link_libraries(cxx_interface INTERFACE -fsanitize=address)
endif()
if(WITH_MSAN)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O1")
  target_compile_options(cxx_interface INTERFACE -fsanitize=memory
                                                 -fno-omit-frame-pointer)
  target_link_libraries(cxx_interface INTERFACE -fsanitize=memory)
endif()
if(WITH_UBSAN)
  target_compile_options(cxx_interface INTERFACE -fsanitize=undefined)
  target_link_libraries(cxx_interface INTERFACE -fsanitize=undefined)
endif()

target_link_libraries(cxx_interface INTERFACE coverage_interface)

#
# Static analysis
#

if(WITH_CLANG_TIDY)
  find_package(ClangTidy "${CMAKE_CXX_COMPILER_VERSION}" EXACT REQUIRED)
  set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=--cuda-host-only")
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
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
  include(FetchContent)
  FetchContent_Declare(
    walberla
    GIT_REPOSITORY https://i10git.cs.fau.de/walberla/walberla
    GIT_TAG        08f04ef64f95609b47838db85862033a1600afa1
    PATCH_COMMAND  patch -p0 -sN --ignore-whitespace -i
                   ${PROJECT_SOURCE_DIR}/cmake/waLBerlaFunctions.patch || echo "Patch not applied!"
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
    add_subdirectory("${walberla_SOURCE_DIR}" "${walberla_BINARY_DIR}")
    set(LB_WALBERLA 1)
    set(WALBERLA_LIBS
        walberla::core walberla::domain_decomposition walberla::blockforest
        walberla::boundary walberla::core walberla::domain_decomposition
        walberla::field walberla::lbm walberla::timeloop walberla::vtk)
    # TODO WALBERLA: remove -lmpi_cxx -lmpi (required to find symbol MPI::Datatype::Free() in VTK code) and -lstdc++fs
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
      set(WALBERLA_LIBS ${WALBERLA_LIBS} stdc++fs -lmpi_cxx -lmpi)
    endif()
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
