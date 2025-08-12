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

cmake_minimum_required(VERSION 3.27.6)
cmake_policy(VERSION 3.27.6)
message(STATUS "CMake version: ${CMAKE_VERSION}")
# CMake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# project info
project(
  ESPResSo
  VERSION "5.0.0"
  LANGUAGES C CXX
  HOMEPAGE_URL "https://espressomd.org"
  DESCRIPTION
    "Extensible Simulation Package for Research on Soft Matter Systems")

# programming language standards
set(ESPRESSO_MINIMAL_C_STANDARD 11)
set(ESPRESSO_MINIMAL_CXX_STANDARD 20)
set(ESPRESSO_MINIMAL_CUDA_STANDARD 20)
set(ESPRESSO_MINIMAL_CUDA_VERSION 12.0)

include(FeatureSummary)
include(GNUInstallDirs)
include(ExternalProject)
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
option(ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM
       "Build with shared memory parallelism support" OFF)
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
# Compiler flags: must be added to all ESPResSo targets
#

add_library(espresso_cpp_flags INTERFACE)
add_library(espresso::cpp_flags ALIAS espresso_cpp_flags)

set_target_properties(
  espresso_cpp_flags PROPERTIES INTERFACE_COMPILE_FEATURES
                                cxx_std_${ESPRESSO_MINIMAL_CXX_STANDARD})

#
# AVX2 support
#

include(CheckCXXCompilerFlag)

add_library(espresso_avx_flags INTERFACE)
add_library(espresso::avx_flags ALIAS espresso_avx_flags)

#
# C/C++ compiler
#

block(SCOPE_FOR VARIABLES)
set(C_COMPILER_IS_OFFICIALLY_TESTED FALSE)
set(CXX_COMPILER_IS_OFFICIALLY_TESTED FALSE)
macro(espresso_minimal_compiler_version)
  foreach(LANG C CXX)
    if(CMAKE_${LANG}_COMPILER_ID STREQUAL "${ARGV0}")
      if(CMAKE_${LANG}_COMPILER_VERSION VERSION_GREATER_EQUAL "${ARGV1}")
        set(${LANG}_COMPILER_IS_OFFICIALLY_TESTED TRUE)
      else()
        message(
          FATAL_ERROR
            "Unsupported compiler ${CMAKE_${LANG}_COMPILER_ID} \
                       ${CMAKE_${LANG}_COMPILER_VERSION} (required version >= ${ARGV1})"
        )
      endif()
    endif()
  endforeach()
endmacro()

espresso_minimal_compiler_version("GNU" 12.2.0)
espresso_minimal_compiler_version("Clang" 18.1.0)
espresso_minimal_compiler_version("AppleClang" 17.0.0)
espresso_minimal_compiler_version("IntelLLVM" 2023.1)

set(UNSUPPORTED_COMPILERS "Intel;MSVC")
foreach(LANG C CXX)
  if(CMAKE_${LANG}_COMPILER_ID IN_LIST UNSUPPORTED_COMPILERS)
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_${LANG}_COMPILER_ID}")
  endif()
  if(NOT ${LANG}_COMPILER_IS_OFFICIALLY_TESTED)
    message(
      AUTHOR_WARNING
        "The ${CMAKE_${LANG}_COMPILER_ID} compiler isn't actively tested by ${PROJECT_NAME}; \
         ${LANG} compiler flags defined by ${PROJECT_NAME}'s CMakeLists.txt might require tweaks"
    )
  endif()
endforeach()
endblock()

#
# CUDA compiler
#

if(ESPRESSO_BUILD_WITH_CUDA)
  include(CheckLanguage)
  enable_language(CUDA)
  check_language(CUDA)
  add_library(espresso_cuda_flags INTERFACE)
  add_library(espresso::cuda_flags ALIAS espresso_cuda_flags)
  find_package(CUDAToolkit ${ESPRESSO_MINIMAL_CUDA_VERSION} REQUIRED)
  if(NOT DEFINED ESPRESSO_CMAKE_CUDA_ARCHITECTURES)
    if("$ENV{CUDAARCHS}" STREQUAL "")
      # 1. sm_61: GTX-1000 series (Pascal)
      # 2. sm_75: RTX-2000 series (Turing)
      # 3. sm_86: RTX-3000 series (Ampere)
      # 4. sm_89: RTX-4000 series (Ada)
      set(ESPRESSO_CUDA_ARCHITECTURES "75;86;89")
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
  set_target_properties(
    espresso_cuda_flags PROPERTIES INTERFACE_COMPILE_FEATURES
                                   cuda_std_${ESPRESSO_MINIMAL_CUDA_STANDARD})
endif()

# Python interpreter and Cython interface library
if(ESPRESSO_BUILD_WITH_PYTHON)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Cython 3.0.8...<3.1.0 REQUIRED)
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
    -Wcast-qual
    -Wcast-align
    -Wunused-macros
    -Wpointer-arith
    -Winit-self
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wformat-signedness>
    $<$<CXX_COMPILER_ID:GNU>:-Wformat=2>
    $<$<CXX_COMPILER_ID:GNU>:-Wduplicated-branches>
    $<$<CXX_COMPILER_ID:GNU>:-Wshadow=compatible-local>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wshadow-field-in-constructor-modified>
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
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:-Wno-format-nonliteral>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:-Wno-float-conversion>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:$<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-implicit-int-float-conversion>>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:$<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-implicit-float-conversion>>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_HDF5}>:-Wno-cast-qual>
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
      # disable warnings from -Wall and -Wextra
      -Wno-array-bounds
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
  unset(SKIP_CLANG_TIDY_CHECKS)
  unset(SKIP_CLANG_TIDY_CHECKS_CXX)
  unset(SKIP_CLANG_TIDY_CHECKS_CUDA)
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

if(ESPRESSO_BUILD_WITH_FFTW)
  if(ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM)
    list(APPEND FFTW3_COMPONENTS omp)
  endif()
  if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
    list(APPEND FFTW3_COMPONENTS mpi)
  endif()
  find_package(fftw3 REQUIRED COMPONENTS ${FFTW3_COMPONENTS})

  if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/heffte-src)
    find_package(Heffte 2.4.1 QUIET)
  endif()
  if(NOT DEFINED Heffte_FOUND OR NOT ${Heffte_FOUND})
    # cmake-format: off
    FetchContent_Declare(
      heffte
      GIT_REPOSITORY https://github.com/icl-utk-edu/heffte.git
      GIT_TAG        v2.4.1
      OVERRIDE_FIND_PACKAGE
    )
    # cmake-format: on
    set(Heffte_ENABLE_FFTW ON CACHE BOOL "")
    if(ESPRESSO_BUILD_WITH_CUDA)
      set(Heffte_ENABLE_CUDA ON CACHE BOOL "")
    endif()
    set(BUILD_SHARED_LIBS ON)
    FetchContent_MakeAvailable(heffte)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    foreach(HEFFTE_MODULE IN ITEMS FFTW CUDA)
      if(TARGET Heffte::${HEFFTE_MODULE})
        set(Heffte_${HEFFTE_MODULE}_FOUND true)
      endif()
    endforeach()
    add_library(Heffte::Heffte INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Heffte::Heffte INTERFACE Heffte)
    install(TARGETS Heffte LIBRARY DESTINATION "${ESPRESSO_INSTALL_LIBDIR}")
    if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
      set(ESPRESSO_P3M_CXX_CLANG_TIDY "${ESPRESSO_CXX_CLANG_TIDY}")
      set(ESPRESSO_P3M_CUDA_CLANG_TIDY "${ESPRESSO_CUDA_CLANG_TIDY}")
      unset(SKIP_CLANG_TIDY_CHECKS)
      unset(SKIP_CLANG_TIDY_CHECKS_CXX)
      unset(SKIP_CLANG_TIDY_CHECKS_CUDA)
      # silence heFFTe diagnostics
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-exception-escape")
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-narrowing-conversions")
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-branch-clone")
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-loop-convert")
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-override")
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
      list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-else-after-return")
      list(APPEND SKIP_CLANG_TIDY_CHECKS
           "-readability-avoid-const-params-in-decls")
      espresso_override_clang_tidy_checks(
        ESPRESSO_P3M_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
        "${SKIP_CLANG_TIDY_CHECKS_CXX}")
      espresso_override_clang_tidy_checks(
        ESPRESSO_P3M_CUDA_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
        "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
    endif()
  endif()
endif()

if(ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM)
  find_package(OpenMP REQUIRED)

  if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/kokkos-src)
    find_package(Kokkos 4.3 QUIET)
  endif()
  if(NOT DEFINED Kokkos_FOUND OR NOT ${Kokkos_FOUND})
    # cmake-format: off
    FetchContent_Declare(
      kokkos
      GIT_REPOSITORY https://github.com/kokkos/kokkos.git
      GIT_TAG        18b830e # version 4.6.1 with patches
      OVERRIDE_FIND_PACKAGE
    )
    # cmake-format: on
    set(BUILD_SHARED_LIBS ON)
    set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
    set(Kokkos_ENABLE_IMPL_VIEW_LEGACY ON CACHE BOOL "")
    set(Kokkos_ENABLE_COMPLEX_ALIGN ON CACHE BOOL "")
    set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON CACHE BOOL "")
    set(Kokkos_ENABLE_HWLOC ON CACHE BOOL "")
    set(Kokkos_ARCH_NATIVE ON CACHE BOOL "")
    FetchContent_MakeAvailable(kokkos)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
    install(TARGETS kokkos LIBRARY DESTINATION "${ESPRESSO_INSTALL_LIBDIR}")
    # install all kokkos shared objects
    get_target_property(ESPRESSO_KOKKOS_LIBS kokkos INTERFACE_LINK_LIBRARIES)
    foreach(target_name IN LISTS ESPRESSO_KOKKOS_LIBS)
      get_target_property(target_type ${target_name} TYPE)
      if(${target_type} STREQUAL "SHARED_LIBRARY" AND ${target_name} MATCHES
                                                      "^kokkos[a-zA-Z0-9_]+$")
        install(TARGETS ${target_name}
                LIBRARY DESTINATION "${ESPRESSO_INSTALL_LIBDIR}")
      endif()
    endforeach()
  endif()

  if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/cabana-src)
    find_package(Cabana 0.7.0 QUIET)
  endif()
  if(NOT DEFINED Cabana_FOUND OR NOT ${Cabana_FOUND})
    # cmake-format: off
    FetchContent_Declare(
      cabana
      GIT_REPOSITORY https://github.com/ECP-copa/Cabana.git
      GIT_TAG        e76c1a1 # 0.7.0 with patches
      PATCH_COMMAND  patch -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cabana.patch
    )
    # cmake-format: on
    set(Cabana_REQUIRE_HEFFTE ${ESPRESSO_BUILD_WITH_FFTW} CACHE BOOL "")
    FetchContent_MakeAvailable(cabana)
  endif()

  if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
    # silence Kokkos and Cabana diagnostics
    set(ESPRESSO_CABANA_CXX_CLANG_TIDY "${ESPRESSO_CXX_CLANG_TIDY}")
    set(ESPRESSO_CABANA_CUDA_CLANG_TIDY "${ESPRESSO_CUDA_CLANG_TIDY}")
    unset(SKIP_CLANG_TIDY_CHECKS)
    unset(SKIP_CLANG_TIDY_CHECKS_CXX)
    unset(SKIP_CLANG_TIDY_CHECKS_CUDA)
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-bugprone-return-const-ref-from-parameter")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-narrowing-conversions")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-sizeof-expression")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-override")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-bool-literals")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-equals-delete")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-equals-default")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-nullptr")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-pass-by-value")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-loop-convert")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-return-braced-init-list")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-else-after-return")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-non-const-parameter")
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-simplify-boolean-expr")
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-readability-avoid-const-params-in-decls")
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-clang-analyzer-optin.performance.Padding")
    list(APPEND SKIP_CLANG_TIDY_CHECKS
         "-clang-analyzer-security.insecureAPI.strcpy")
    espresso_override_clang_tidy_checks(
      ESPRESSO_CABANA_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
      "${SKIP_CLANG_TIDY_CHECKS_CXX}")
    espresso_override_clang_tidy_checks(
      ESPRESSO_CABANA_CUDA_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
      "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
    if(ESPRESSO_BUILD_WITH_FFTW)
      espresso_override_clang_tidy_checks(
        ESPRESSO_P3M_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
        "${SKIP_CLANG_TIDY_CHECKS_CXX}")
      espresso_override_clang_tidy_checks(
        ESPRESSO_P3M_CUDA_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
        "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
    endif()
  endif()
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
  add_library(hdf5 INTERFACE)
  target_link_libraries(hdf5 INTERFACE $<BUILD_INTERFACE:${HDF5_LIBRARIES}>)
  target_include_directories(hdf5
                             INTERFACE $<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>)
  # cmake-format: off
  FetchContent_Declare(
    HighFive
    GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
    GIT_TAG        0103467 # v3.0.0-beta with patches
  )
  # cmake-format: on
  FetchContent_MakeAvailable(HighFive)
endif()

if(ESPRESSO_BUILD_WITH_SCAFACOS)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SCAFACOS scafacos REQUIRED)
endif()

if(ESPRESSO_BUILD_WITH_GSL)
  find_package(GSL REQUIRED)
endif()

if(ESPRESSO_BUILD_WITH_STOKESIAN_DYNAMICS)
  # cmake-format: off
  FetchContent_Declare(
    stokesian_dynamics
    GIT_REPOSITORY https://github.com/hmenke/espresso-stokesian-dynamics.git
    GIT_TAG        90f6de70d1c0ac9cf468f2fc90f9c601139e4c88
  )
  # cmake-format: on
  set(STOKESIAN_DYNAMICS 1)
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
  FetchContent_MakeAvailable(stokesian_dynamics)
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_OLD_CMAKE_INSTALL_LIBDIR}")
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

if(POLICY CMP0167)
  # use BoostConfig.cmake shipped with Boost 1.70+ instead of the one in CMake
  cmake_policy(SET CMP0167 NEW)
endif()

list(APPEND ESPRESSO_BOOST_COMPONENTS mpi serialization)

if(ESPRESSO_BUILD_TESTS)
  list(APPEND ESPRESSO_BOOST_COMPONENTS unit_test_framework)
endif()

find_package(Boost 1.83.0 REQUIRED ${ESPRESSO_BOOST_COMPONENTS})

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
  add_library(espresso_tests_cpp_flags INTERFACE)
  add_library(espresso::tests::cpp_flags ALIAS espresso_tests_cpp_flags)
  target_compile_options(espresso_tests_cpp_flags INTERFACE -Wno-unused-macros)
  if(ESPRESSO_BUILD_WITH_COVERAGE)
    target_compile_options(
      espresso_tests_cpp_flags
      INTERFACE $<$<CXX_COMPILER_ID:GNU>:-fno-default-inline>
                $<$<CXX_COMPILER_ID:GNU>:-fno-elide-constructors>)
  endif()
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
  # cmake-format: off
  FetchContent_Declare(
    walberla
    GIT_REPOSITORY https://i10git.cs.fau.de/walberla/walberla.git
    GIT_TAG        c69cb11d # v7.1 with patches
  )
  # cmake-format: on
  string(REGEX REPLACE "([/\\]walberla)-src$" "\\1-build" walberla_BINARY_DIR
                       "${walberla_SOURCE_DIR}")
  set(WALBERLA_BUILD_TESTS off CACHE BOOL "")
  set(WALBERLA_BUILD_BENCHMARKS off CACHE BOOL "")
  set(WALBERLA_BUILD_TOOLS off CACHE BOOL "")
  set(WALBERLA_BUILD_TUTORIALS off CACHE BOOL "")
  set(WALBERLA_BUILD_SHOWCASES off CACHE BOOL "")
  set(WALBERLA_BUILD_DOC off CACHE BOOL "")
  set(WALBERLA_LOGLEVEL "WARNING" CACHE STRING "")
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
  if(ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM)
    set(WALBERLA_BUILD_WITH_OPENMP on CACHE BOOL "")
  endif()
  set(WALBERLA_BUILD_WITH_FASTMATH off CACHE BOOL "")
  set(BUILD_SHARED_LIBS ON)
  set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
  FetchContent_MakeAvailable(walberla)
  set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
  set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
  set(WALBERLA_LIBS
      walberla::core walberla::domain_decomposition walberla::blockforest
      walberla::communication walberla::field walberla::stencil walberla::vtk)
  if(WALBERLA_BUILD_WITH_FFTW)
    list(APPEND WALBERLA_LIBS walberla::fft)
  endif()
  if(WALBERLA_BUILD_WITH_CUDA)
    list(APPEND WALBERLA_LIBS walberla::gpu)
  endif()
  set(WALBERLA_LIBS_EXPORT
      ${WALBERLA_LIBS} lodepng walberla::boundary walberla::executiontree
      walberla::geometry walberla::lbm walberla::timeloop)
  # install all waLBerla shared objects
  foreach(target_alias IN LISTS WALBERLA_LIBS_EXPORT)
    string(REPLACE "::" "_" target_name ${target_alias})
    get_target_property(target_type ${target_name} TYPE)
    if(${target_type} STREQUAL "SHARED_LIBRARY")
      install(TARGETS ${target_name}
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
  # cmake-format: off
  FetchContent_Declare(
    caliper
    GIT_REPOSITORY https://github.com/LLNL/Caliper.git
    GIT_TAG        v2.13.0
  )
  # cmake-format: on
  set(CALIPER_OPTION_PREFIX on CACHE BOOL "")
  set(CALIPER_WITH_MPI on CACHE BOOL "")
  set(CALIPER_WITH_NVTX off CACHE BOOL "")
  set(CALIPER_WITH_CUPTI off CACHE BOOL "")
  set(CALIPER_INSTALL_CONFIG off CACHE BOOL "")
  set(CALIPER_INSTALL_HEADERS off CACHE BOOL "")
  set(BUILD_SHARED_LIBS ON)
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
  set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
  FetchContent_MakeAvailable(caliper)
  set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_OLD_CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
endif()

#
# Set source file properties
#

function(espresso_set_common_target_properties)
  set(TARGET_NAME "${ARGV0}")
  get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)
  if(NOT TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
    # set non-transitive properties; transitive properties should be propagated
    # to dependent targets via INTERFACE targets, such as espresso::cpp_flags
    set_target_properties(${TARGET_NAME} PROPERTIES C_EXTENSIONS OFF)
    set_target_properties(${TARGET_NAME} PROPERTIES CXX_EXTENSIONS OFF)
    set_target_properties(${TARGET_NAME} PROPERTIES CUDA_EXTENSIONS OFF)
    set_target_properties(
      ${TARGET_NAME} PROPERTIES CXX_CLANG_TIDY "${ESPRESSO_CXX_CLANG_TIDY}")
    set_target_properties(
      ${TARGET_NAME} PROPERTIES CUDA_CLANG_TIDY "${ESPRESSO_CUDA_CLANG_TIDY}")
  endif()
endfunction()

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
