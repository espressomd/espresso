#
# Copyright (C) 2009-2025 The ESPResSo project
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
if(POLICY CMP0169)
  # use deprecated FetchContent_Populate function for h5xx dependency
  cmake_policy(SET CMP0169 OLD)
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

espresso_minimal_compiler_version("GNU" 12.2.0)
espresso_minimal_compiler_version("Clang" 18.1.0)
espresso_minimal_compiler_version("AppleClang" 16.0.0)
espresso_minimal_compiler_version("IntelLLVM" 2023.1)

set(UNSUPPORTED_COMPILERS "Intel;MSVC")
if(CMAKE_CXX_COMPILER_ID IN_LIST UNSUPPORTED_COMPILERS)
  message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

include(FeatureSummary)
include(GNUInstallDirs)
include(FetchContent)
include(espresso_option_enum)
include(espresso_enable_avx2_support)
include(espresso_override_clang_tidy_checks)

if(EXISTS "${PROJECT_BINARY_DIR}/CMakeLists.txt")
  message(
    FATAL_ERROR
      "${PROJECT_NAME} cannot be built in-place. Instead, create a build directory and run CMake from there. A new file 'CMakeCache.txt' and a new folder 'CMakeFiles' have just been created by CMake in the current folder and need to be removed."
  )
endif()

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

# build targets as static libraries unless otherwise specified
set(ESPRESSO_BUILD_SHARED_LIBS_DEFAULT OFF)
set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})

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
  GIT_TAG        f36fa0a68bae59f0b516f6587ea8fa7c24a41141
  PATCH_COMMAND  patch -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/cmake/walberla.patch
)
FetchContent_Declare(
  stokesian_dynamics
  GIT_REPOSITORY https://github.com/hmenke/espresso-stokesian-dynamics.git
  GIT_TAG        862a7537a366f0c32f0c25e46bd107bea590faea
)
FetchContent_Declare(
  h5xx
  GIT_REPOSITORY https://github.com/h5md/h5xx.git
  GIT_TAG        0.9.1
)
FetchContent_Declare(
  caliper
  GIT_REPOSITORY https://github.com/LLNL/Caliper.git
  GIT_TAG        v2.12.1
)
# cmake-format: on

# ##############################################################################
# User input options
# ##############################################################################

option(ESPRESSO_BUILD_WITH_CABANA "Build with Cabana" ON)

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
option(ESPRESSO_BUILD_WITH_FPE
       "Build with floating-point exceptions instrumentation" OFF)
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
    find_package(CUDACompilerClang 18.1.0 REQUIRED)
  else()
    message(FATAL_ERROR "Unknown CUDA compiler '${CMAKE_CUDA_COMPILER_ID}'")
  endif()
endif()

# Python interpreter and Cython interface library
if(ESPRESSO_BUILD_WITH_PYTHON)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Cython 0.29.28...<3.0.12 REQUIRED)
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

set(ESPRESSO_OLD_CMAKE_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")

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
# Compiler flags: must be added to all ESPResSo targets
#

add_library(espresso_cpp_flags INTERFACE)
add_library(espresso::cpp_flags ALIAS espresso_cpp_flags)

#
# Compiler diagnostics
#

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
    $<$<CXX_COMPILER_ID:Clang>:-Wnon-c-typedef-for-linkage>
    -Wdelete-non-virtual-dtor
    # disable warnings from -Wall and -Wextra
    -Wno-sign-compare
    -Wno-unused-function
    -Wno-unused-parameter
    -Wno-array-bounds
    $<$<CXX_COMPILER_ID:GNU>:-Wno-restrict>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-clobbered>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-cast-function-type>
    $<$<CXX_COMPILER_ID:IntelLLVM>:-diag-disable=592>
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-gnu-zero-variadic-macro-arguments>
    $<$<NOT:$<CXX_COMPILER_ID:IntelLLVM>>:-Wno-implicit-fallthrough>
    $<$<NOT:$<CXX_COMPILER_ID:IntelLLVM,GNU>>:-Wno-unused-private-field>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm>>:-Wno-psabi>
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
      $<$<NOT:$<CXX_COMPILER_ID:IntelLLVM>>:-Wno-implicit-fallthrough>
      # warnings are errors
      $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<CUDA_COMPILER_ID:NVIDIA>>:--Werror=all-warnings>
      $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<CUDA_COMPILER_ID:Clang>>:-Werror>
  )
endif()

#
# Code coverage
#

if(ESPRESSO_BUILD_WITH_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
      espresso_cpp_flags INTERFACE -g -fprofile-instr-generate
                                   -fcoverage-mapping)
  else()
    target_compile_options(espresso_cpp_flags INTERFACE -g --coverage
                                                        -fprofile-abs-path)
    target_link_libraries(espresso_cpp_flags INTERFACE gcov)
  endif()
endif()

#
# Portability options
#

# prevent 80-bit arithmetic on old Intel processors
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SIZEOF_VOID_P EQUAL 4
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "[xX]86")
  target_compile_options(espresso_cpp_flags INTERFACE -ffloat-store)
endif()

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

#
# Static analysis
#

if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
  find_package(ClangTidy "${CMAKE_CXX_COMPILER_VERSION}" EXACT REQUIRED)
  set(ESPRESSO_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  set(ESPRESSO_CUDA_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=--cuda-host-only")
  set(SKIP_CLANG_TIDY_CHECKS "")
  set(SKIP_CLANG_TIDY_CHECKS_CXX "")
  set(SKIP_CLANG_TIDY_CHECKS_CUDA "")
  if(ESPRESSO_BUILD_WITH_CALIPER)
    # Clang-Tidy sometimes emits diagnostics in code enclosed in `extern "C"`
    # that are not always actionable, since they may rely on keywords only
    # available in the C++ language. While some checks have an extra flag
    # 'IgnoreExternC' to disable them inside C code, not all affected checks
    # have been fixed yet. For an in-depth discussion on this topic, see
    # https://github.com/llvm/llvm-project/issues/35272
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-nullptr")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-deprecated-headers")
  endif()
  if(ESPRESSO_BUILD_WITH_CUDA)
    # silence casts in cuda_runtime.h (for both C++ and CUDA source files)
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-casting-through-void")
    # silence nullptr dereference in cuda::thrust
    list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA
         "-clang-analyzer-core.NonNullParamChecker")
  endif()
  espresso_override_clang_tidy_checks(
    ESPRESSO_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_CXX}")
  espresso_override_clang_tidy_checks(
    ESPRESSO_CUDA_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
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
# Libraries
#

if(ESPRESSO_BUILD_WITH_CABANA)
  # TODO Kokkos_DIR needs to be set here

  find_package(Kokkos REQUIRED)

  include(FetchContent)
  if(EXISTS "${CMAKE_BINARY_DIR}/_deps/cabana-build")
    FetchContent_Declare(cabana
      GIT_REPOSITORY https://github.com/ECP-copa/Cabana.git
      GIT_TAG 07573782590e223edc6e5f83ec0807be5941df14
    )
  else()# Fix cabana doxygen conflict until it is fixed upstream
    FetchContent_Declare(cabana
      GIT_REPOSITORY https://github.com/ECP-copa/Cabana.git
      GIT_TAG 07573782590e223edc6e5f83ec0807be5941df14
      PATCH_COMMAND patch -N --reject-file=- -i "${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_cabana_doxygen.diff" 
    )
  endif()
  FetchContent_MakeAvailable(cabana)
endif()

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
    if(${HDF5_VERSION} VERSION_GREATER_EQUAL 1.13.0)
      message(
        WARNING
          "Library hdf5 >= 1.13.0 is incompatible with h5xx 0.9.1 on some platforms (https://github.com/fhoefling/h5xx/issues/11)"
      )
    endif()
  endif()
  # We need to define our own h5xx target, since the project version we depend
  # on is not fully compatible with modern CMake.
  FetchContent_GetProperties(h5xx)
  if(NOT h5xx_POPULATED)
    FetchContent_Populate(h5xx)
  endif()
  add_library(h5xx INTERFACE)
  add_library(hdf5 INTERFACE)
  target_link_libraries(hdf5 INTERFACE $<BUILD_INTERFACE:${HDF5_LIBRARIES}>)
  target_include_directories(hdf5
                             INTERFACE $<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>)
  target_include_directories(h5xx
                             INTERFACE $<BUILD_INTERFACE:${h5xx_SOURCE_DIR}>)
  target_compile_features(h5xx INTERFACE cxx_std_11)
  target_compile_definitions(h5xx INTERFACE H5XX_USE_MPI)
  add_library(espresso_h5xx_cpp_flags INTERFACE)
  add_library(espresso::h5xx_cpp_flags ALIAS espresso_h5xx_cpp_flags)
  target_compile_options(
    espresso_h5xx_cpp_flags
    INTERFACE
      $<$<CXX_COMPILER_ID:GNU>:-Wno-terminate>
      $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-exceptions>
      $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-unused-but-set-variable>
      $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-delete-non-abstract-non-virtual-dtor>
  )
  if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
    set(H5XX_CXX_CLANG_TIDY "${ESPRESSO_CXX_CLANG_TIDY}")
    set(SKIP_CLANG_TIDY_CHECKS "")
    set(SKIP_CLANG_TIDY_CHECKS_CXX "")
    # silence hdf5 and h5xx diagnostics
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-clang-analyzer-deadcode.DeadStores")
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-clang-analyzer-optin.performance.Padding")
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-clang-diagnostic-unused-but-set-variable")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-clang-diagnostic-exceptions")
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-bugprone-multi-level-implicit-pointer-conversion")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-assignment-in-if-condition")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-exception-escape")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-narrowing-conversions")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-nullptr")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-noexcept")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-override")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-equals-default")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-loop-convert")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-pass-by-value")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-redundant-void-arg")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-else-after-return")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-container-size-empty")
    espresso_override_clang_tidy_checks(
      H5XX_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
      "${SKIP_CLANG_TIDY_CHECKS_CXX}")
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
  FetchContent_GetProperties(stokesian_dynamics)
  set(STOKESIAN_DYNAMICS 1)
  if(NOT stokesian_dynamics_POPULATED)
    set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
    FetchContent_MakeAvailable(stokesian_dynamics)
    set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_OLD_CMAKE_INSTALL_LIBDIR}")
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
set(ESPRESSO_MPIEXEC_PREFLAGS "")
# Open MPI 4.x has a bug on NUMA archs that prevents running in singleton mode
set(ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA OFF)
set(ESPRESSO_CPU_MODEL_NAME_OMPI_SINGLETON_NUMA_PATTERN "AMD (EPYC|Ryzen)")

if("${ESPRESSO_MPIEXEC_VENDOR}" STREQUAL "OpenMPI")
  list(APPEND ESPRESSO_MPIEXEC_PREFLAGS "--oversubscribe")
  if(ESPRESSO_INSIDE_DOCKER)
    list(APPEND ESPRESSO_MPIEXEC_PREFLAGS "--bind-to" "none")
  endif()
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

#
# Boost
#

list(APPEND BOOST_COMPONENTS mpi serialization filesystem system)

if(ESPRESSO_BUILD_TESTS)
  list(APPEND BOOST_COMPONENTS unit_test_framework)
endif()

find_package(Boost 1.74.0 REQUIRED ${BOOST_COMPONENTS})

# Patches for older Boost releases compiled with C++17 or above; while
# originally meant for MSVC, the Clang compiler is also affected.
# https://boost-users.boost.narkive.com/h2eXAUzU
# https://devblogs.microsoft.com/cppblog/c17-feature-removals-and-deprecations/
# https://github.com/PixarAnimationStudios/OpenUSD/issues/2634
target_compile_options(
  espresso_cpp_flags
  INTERFACE
    $<$<VERSION_LESS:${Boost_VERSION},1.83.0>:-D_HAS_AUTO_PTR_ETC=0>
    $<$<VERSION_LESS:${Boost_VERSION},1.81.0>:-DBOOST_NO_CXX98_FUNCTION_BASE>)

# enable boost::variant with more than 20 types
target_compile_options(
  espresso_cpp_flags INTERFACE -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
                               -DBOOST_MPL_LIMIT_LIST_SIZE=30)

#
# Paths
#

set(CMAKE_INSTALL_RPATH "${ESPRESSO_INSTALL_LIBDIR}")

#
# Packaging
#

# drop 'lib' prefix from all libraries
set(ESPRESSO_SHARED_LIBRARY_PREFIX "")
set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")

set(CMAKE_MACOSX_RPATH TRUE)

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
    set(BUILD_SHARED_LIBS ON)
    set(CMAKE_SHARED_LIBRARY_PREFIX "libwalberla_")
    FetchContent_MakeAvailable(walberla)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
  endif()
  set(WALBERLA_LIBS
      walberla::core walberla::domain_decomposition walberla::blockforest
      walberla::boundary walberla::field walberla::lbm walberla::timeloop
      walberla::vtk)
  if(WALBERLA_BUILD_WITH_FFTW)
    list(APPEND WALBERLA_LIBS walberla::fft)
  endif()
  if(WALBERLA_BUILD_WITH_CUDA)
    list(APPEND WALBERLA_LIBS walberla::gpu)
  endif()
  set(WALBERLA_LIBS_EXPORT
      ${WALBERLA_LIBS}
      walberla::lodepng
      walberla::sqlite3
      walberla::lbm_generated
      walberla::pe
      walberla::pe_coupling
      walberla::pde
      walberla::mesa_pd
      walberla::lbm_mesapd_coupling
      walberla::sqlite
      walberla::executiontree
      walberla::geometry
      walberla::postprocessing
      walberla::gather
      walberla::gui)
  # workaround for https://gitlab.kitware.com/cmake/cmake/-/issues/21283
  foreach(target_w_namespace IN LISTS WALBERLA_LIBS)
    string(REPLACE "walberla::" "" target_wo_namespace ${target_w_namespace})
    add_library(${target_w_namespace} ALIAS ${target_wo_namespace})
  endforeach()
  # install all waLBerla shared objects
  foreach(target_w_namespace IN LISTS WALBERLA_LIBS_EXPORT)
    string(REPLACE "walberla::" "" target_wo_namespace ${target_w_namespace})
    get_target_property(target_type ${target_wo_namespace} TYPE)
    if(${target_type} STREQUAL "SHARED_LIBRARY")
      install(TARGETS ${target_wo_namespace}
              LIBRARY DESTINATION "${ESPRESSO_INSTALL_LIBDIR}")
    endif()
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
  set(CALIPER_INSTALL_CONFIG off CACHE BOOL "")
  set(CALIPER_INSTALL_HEADERS off CACHE BOOL "")
  if(NOT caliper_POPULATED)
    set(BUILD_SHARED_LIBS ON)
    set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
    set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
    FetchContent_MakeAvailable(caliper)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_OLD_CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
  endif()
  target_compile_options(caliper-common
                         PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-restrict>)
  target_compile_options(
    caliper-services
    PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wno-deprecated-declarations>)
  target_compile_options(
    caliper-runtime
    PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
            $<$<CXX_COMPILER_ID:GNU>:-Wno-volatile>
            $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-deprecated-volatile>)
  set_target_properties(caliper-runtime PROPERTIES CXX_CLANG_TIDY "")
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
