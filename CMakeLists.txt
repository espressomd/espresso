#
# Copyright (C) 2009-2026 The ESPResSo project
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
if(ESPRESSO_BUILD_WITH_CUDA)
  cmake_minimum_required(VERSION 3.31.0)
endif()
cmake_policy(VERSION 3.27.6)
if(POLICY CMP0167)
  # use BoostConfig.cmake shipped with Boost 1.70+ instead of the one in CMake
  cmake_policy(SET CMP0167 NEW)
endif()
message(STATUS "CMake version: ${CMAKE_VERSION}")
# CMake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# project info
project(
  ESPResSo
  VERSION "5.0.1"
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

# shared objects require position-independent code in static libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# On Mac OS X, first look for other packages, then frameworks
set(CMAKE_FIND_FRAMEWORK LAST)

# avoid applying patches twice with FetchContent
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

# customize default options based on basic information
set(ESPRESSO_BUILD_WITH_WALBERLA_AVX_DEFAULT OFF)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "[xX]86")
  set(ESPRESSO_BUILD_WITH_WALBERLA_AVX_DEFAULT ON)
endif()

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
option(ESPRESSO_BUILD_WITH_NLOPT "Build with NLopt support" OFF)
option(ESPRESSO_BUILD_WITH_WALBERLA "Build with waLBerla support" ON)
option(ESPRESSO_BUILD_WITH_WALBERLA_AVX
       "Build waLBerla kernels with AVX2 vectorization"
       ${ESPRESSO_BUILD_WITH_WALBERLA_AVX_DEFAULT})
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
if(ESPRESSO_BUILD_WITH_PYTHON)
  set(ESPRESSO_MODULE_INSTALL_PATH ""
      CACHE PATH "Install ESPResSo Python module in a custom location")
endif()

if(ESPRESSO_BUILD_WITH_CCACHE)
  find_program(CCACHE_PROGRAM ccache REQUIRED)
  if(CCACHE_PROGRAM)
    message(STATUS "Found ccache: ${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_CUDA_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
  endif()
endif()

set(ESPRESSO_NUMPY_SITEARCH ""
    CACHE FILEPATH
          "NumPy's third-party platform dependent installation directory")

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

add_library(espresso_compiler_flags INTERFACE)
add_library(espresso::compiler_flags ALIAS espresso_compiler_flags)

set_property(
  TARGET espresso_compiler_flags APPEND
  PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_${ESPRESSO_MINIMAL_CXX_STANDARD})

add_library(espresso_kokkos_compiler_flags INTERFACE)
add_library(espresso::kokkos::compiler_flags ALIAS
            espresso_kokkos_compiler_flags)

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
espresso_minimal_compiler_version("CrayClang" 17.0.0)
espresso_minimal_compiler_version("IntelLLVM" 2023.1)
espresso_minimal_compiler_version("NVHPC" 25.5)

set(ESPRESSO_UNSUPPORTED_COMPILERS "Intel;MSVC")
foreach(LANG C CXX)
  if(CMAKE_${LANG}_COMPILER_ID IN_LIST ESPRESSO_UNSUPPORTED_COMPILERS)
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

set(ESPRESSO_XCOMPILER "")
set(ESPRESSO_PTXAS "")
if(ESPRESSO_BUILD_WITH_CUDA)
  include(CheckLanguage)
  enable_language(CUDA)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVHPC")
    find_package(CUDAToolkit REQUIRED)
  else()
    find_package(CUDAToolkit ${ESPRESSO_MINIMAL_CUDA_VERSION} REQUIRED)
  endif()
  if(NOT DEFINED ESPRESSO_CMAKE_CUDA_ARCHITECTURES)
    if("$ENV{CUDAARCHS}" STREQUAL "")
      # cmake-format: off
      # sm_60: P100 series (Pascal)
      # sm_61: GTX-1000 series (Pascal)
      # sm_70: V100 series (Volta)
      # sm_75: RTX-2000 series (Turing)
      # sm_86: RTX-3000 series (Ampere)
      # sm_89: RTX-4000 series (Ada)
      # sm_90: H100 series (Hopper)
      # sm_120: RTX-5000 series (Blackwell)
      # cmake-format: on
      set(ESPRESSO_CUDA_ARCHITECTURES "75;86;89")
      if(Kokkos_ENABLE_CUDA)
        set(ESPRESSO_CUDA_ARCHITECTURES "86")
      endif()
    else()
      set(ESPRESSO_CUDA_ARCHITECTURES "$ENV{CUDAARCHS}")
    endif()
    set(ESPRESSO_CMAKE_CUDA_ARCHITECTURES "${ESPRESSO_CUDA_ARCHITECTURES}"
        CACHE INTERNAL "")
  endif()
  set(CMAKE_CUDA_ARCHITECTURES "${ESPRESSO_CMAKE_CUDA_ARCHITECTURES}")
  list(LENGTH CMAKE_CUDA_ARCHITECTURES ESPRESSO_CMAKE_CUDA_ARCHITECTURES_LEN)
  if(Kokkos_ENABLE_CUDA AND ESPRESSO_CMAKE_CUDA_ARCHITECTURES_LEN GREATER 1
     AND NOT CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
    message(
      FATAL_ERROR
        "Kokkos can only be built against a single CUDA architecture when using the nvcc_wrapper, got CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}"
    )
  endif()
  cmake_path(GET CUDA_cuda_driver_LIBRARY PARENT_PATH ESPRESSO_LIBCUDA_RPATH)
  cmake_path(GET CUDA_cudart_LIBRARY PARENT_PATH ESPRESSO_LIBCUDART_RPATH)
  macro(espresso_add_cuda_rpaths)
    set(TARGET_NAME ${ARGV0})
    set_property(
      TARGET ${TARGET_NAME} APPEND
      PROPERTY BUILD_RPATH "${ESPRESSO_LIBCUDA_RPATH}"
               "${ESPRESSO_LIBCUDART_RPATH}")
    set_property(
      TARGET ${TARGET_NAME} APPEND
      PROPERTY INSTALL_RPATH "${ESPRESSO_LIBCUDA_RPATH}"
               "${ESPRESSO_LIBCUDART_RPATH}")
  endmacro()
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    find_package(CUDACompilerNVCC ${ESPRESSO_MINIMAL_CUDA_VERSION} REQUIRED)
    set(ESPRESSO_XCOMPILER "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=>")
    set(ESPRESSO_PTXAS "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xptxas=>")
  elseif(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
    if(ESPRESSO_BUILD_WITH_COVERAGE)
      message(
        FATAL_ERROR
          "Cannot enable code coverage with Clang as the CUDA compiler")
    endif()
    find_package(CUDACompilerClang 18.1.0 REQUIRED)
  elseif(CMAKE_CUDA_COMPILER_ID STREQUAL "NVHPC")
    if(ESPRESSO_BUILD_WITH_COVERAGE)
      message(
        FATAL_ERROR
          "Cannot enable code coverage with nvc++ as the CUDA compiler")
    endif()
    find_package(CUDACompilerNVHPC 25.5 REQUIRED)
    if(${ESPRESSO_NVHPC_DETECTED_CUDA_VERSION} VERSION_LESS
       ${ESPRESSO_MINIMAL_CUDA_VERSION})
      message(
        FATAL_ERROR
          "Minimal CUDA version is ${ESPRESSO_MINIMAL_CUDA_VERSION}, but nvc++ "
          "reported ${ESPRESSO_NVHPC_DETECTED_CUDA_VERSION}; if your HPC SDK "
          "ships multiple CUDA toolkits, try hinting the most recent one with "
          "'-D CUDAToolkit_NVCC_EXECUTABLE=\${NVHPC_ROOT}/compilers/bin/nvcc'")
    endif()
  else()
    message(FATAL_ERROR "Unknown CUDA compiler '${CMAKE_CUDA_COMPILER_ID}'")
  endif()
  set_property(
    TARGET espresso_compiler_flags APPEND
    PROPERTY INTERFACE_COMPILE_FEATURES
             cuda_std_${ESPRESSO_MINIMAL_CUDA_STANDARD})
endif()

# Python interpreter and Cython interface library
if(ESPRESSO_BUILD_WITH_PYTHON)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Cython 3.0.4...<4.0.0 REQUIRED)
  if(CYTHON_VERSION VERSION_LESS 3.0.8)
    message(WARNING "We strongly recommend using Cython 3.0.8 or later")
  endif()
  find_program(IPYTHON_EXECUTABLE NAMES jupyter ipython3 ipython)
  if(NOT DEFINED CACHE{ESPRESSO_NUMPY_SITEARCH}
     OR "$CACHE{ESPRESSO_NUMPY_SITEARCH}" STREQUAL "")
    # NumPy isn't necessarily installed in Python_SITEARCH. EasyBuild packages
    # NumPy within SciPy-bundle, for example. Cython needs to include NumPy's
    # parent directory to properly detect file numpy/__init__.pxd.
    block(SCOPE_FOR VARIABLES PROPAGATE ESPRESSO_NUMPY_SITEARCH)
    execute_process(
      COMMAND ${Python_EXECUTABLE} -c "import numpy;print(numpy.__file__)"
      OUTPUT_VARIABLE numpy_init_path ERROR_VARIABLE numpy_init_error
      RESULT_VARIABLE numpy_init_result)
    if(NOT ${numpy_init_result} EQUAL 0)
      message(FATAL_ERROR "Cannot import numpy:\n${numpy_init_error}")
    endif()
    cmake_path(GET numpy_init_path PARENT_PATH numpy_path)
    cmake_path(GET numpy_path PARENT_PATH numpy_sitearch)
    set(ESPRESSO_NUMPY_SITEARCH "${numpy_sitearch}"
        CACHE FILEPATH
              "NumPy's third-party platform dependent installation directory"
              FORCE)
    endblock()
  endif()
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

# python site-packages, can be overridden with CMake options
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
  if(NOT ESPRESSO_MODULE_INSTALL_PATH STREQUAL "")
    set(ESPRESSO_INSTALL_PYTHON "${ESPRESSO_MODULE_INSTALL_PATH}")
  endif()
  # install shared objects in the Python folder
  set(ESPRESSO_INSTALL_LIBDIR "${ESPRESSO_INSTALL_PYTHON}/espressomd")
  add_custom_target(espresso_packaging_dependencies)
endif()

#
# Compiler diagnostics
#

# cmake-format: off
target_compile_options(
  espresso_compiler_flags
  INTERFACE
    $<$<COMPILE_LANGUAGE:CXX,CUDA>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX,CUDA>:-Wextra>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang>:-pedantic>
    # add extra warnings
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-Wfloat-conversion>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-Wdelete-non-virtual-dtor>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-Wnon-virtual-dtor>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-Wunused-macros>
    $<$<COMPILE_LANGUAGE:CXX>:-Wcast-qual>
    $<$<COMPILE_LANGUAGE:CXX>:-Wcast-align>
    $<$<COMPILE_LANGUAGE:CXX>:-Wpointer-arith>
    $<$<COMPILE_LANGUAGE:CXX>:-Winit-self>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-Wextern-initializer>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-Wrange-loop-analysis>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-Wnon-c-typedef-for-linkage>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-Wshadow-uncaptured-local>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,IntelLLVM>:-Wimplicit-float-conversion>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,IntelLLVM>:-Wunused-exception-parameter>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,IntelLLVM>:-Wmissing-variable-declarations>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wctad-maybe-unsupported>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wbuiltin-macro-redefined>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wformat-signedness>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wdiv-by-zero>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wextra-semi>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wexceptions>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wdeprecated-declarations>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wdeprecated-enum-enum-conversion>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wdeprecated-enum-float-conversion>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wformat=2>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wbidi-chars>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wcomma-subscript>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wduplicated-branches>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wshadow=compatible-local>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,IntelLLVM>:-Wshadow-field-in-constructor-modified>
    # disable warnings from -Wall and -Wextra, as well as warnings triggered by third-party libraries
    $<$<COMPILE_LANGUAGE:CXX,CUDA>:-Wno-sign-compare>
    $<$<COMPILE_LANGUAGE:CXX,CUDA>:-Wno-unused-parameter>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-Wno-array-bounds>
    $<$<COMPILE_LANG_AND_ID:CUDA,Clang>:-Wno-array-bounds>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<IN_LIST:${CUDA_HOST_COMPILER_ID},GNU;Clang;CrayClang;AppleClang;IntelLLVM>>:-Xcompiler=-Wno-array-bounds>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-restrict>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-cast-function-type>
    $<$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>:-diag-disable=592>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,IntelLLVM>:-Wno-global-constructors>
    $<$<AND:$<BOOL:${ESPRESSO_BUILD_WITH_FFTW}>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>>:-Wno-extra-semi>
    $<$<AND:$<BOOL:${ESPRESSO_BUILD_WITH_HDF5}>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>>:-Wno-extra-semi>
    $<$<AND:$<BOOL:${ESPRESSO_BUILD_WITH_HDF5}>,$<COMPILE_LANGUAGE:CXX>>:-Wno-cast-qual>
    $<$<AND:$<BOOL:${ESPRESSO_BUILD_WITH_HDF5}>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>>:-Wno-old-style-cast>
    $<$<AND:$<BOOL:${ESPRESSO_BUILD_WITH_HDF5}>,$<COMPILE_LANG_AND_ID:CXX,NVHPC>>:--diag_suppress=storage_class_not_first>
    $<$<AND:$<VERSION_LESS:${ESPRESSO_MINIMAL_CXX_STANDARD},23>,$<COMPILE_LANG_AND_ID:CXX,NVHPC>>:--diag_suppress=bad_pp_directive_keyword>
    $<$<COMPILE_LANG_AND_ID:CXX,NVHPC,NVIDIA>:--diag_suppress=code_is_unreachable>
    $<$<COMPILE_LANG_AND_ID:CXX,NVHPC,NVIDIA>:--diag_suppress=loop_not_reachable>
    # disable GCC's C++17 processor-specific ABI warning (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94383)
    $<$<AND:$<COMPILE_LANG_AND_ID:CXX,GNU>,$<IN_LIST:${CMAKE_SYSTEM_PROCESSOR},arm;arm64;aarch64;aarch64_be;powerpc64le>>:-Wno-psabi>
    # warnings are errors
    $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<COMPILE_LANGUAGE:CXX>>:-Werror>
    $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:--Werror=all-warnings>
    $<$<AND:$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>,$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>>:-Werror>
)
# cmake-format: on

#
# Compiler optimizations
#

# cmake-format: off
target_compile_options(
  espresso_compiler_flags
  INTERFACE
    # prevent loop auto-vectorization when doing so can introduce invalid math;
    # for example, when a loop body contains a conditional to skip divisions by
    # zero, the conditional can be ignored (speculative execution) to get better
    # performance with SIMD types, and the invalid values are masked afterwards
    # (but the SIGFPE signal is still sent as a side effect)
    $<$<AND:$<BOOL:${ESPRESSO_BUILD_WITH_FPE}>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>>:-ftrapping-math>
    # configurations for NVCC
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CONFIG:Debug>>:-g -G>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CONFIG:Release>>:-Xptxas=-O3 -Xcompiler=-O3 -DNDEBUG>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CONFIG:MinSizeRel>>:-Xptxas=-O2 -Xcompiler=-Os -DNDEBUG>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CONFIG:RelWithDebInfo>>:-Xptxas=-O2 -Xcompiler=-O2,-g -DNDEBUG>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CONFIG:Coverage>>:-Xptxas=-O3 -Xcompiler=-Og,-g>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CONFIG:RelWithAssert>>:-Xptxas=-O3 -Xcompiler=-O3,-g>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<BOOL:${CMAKE_OSX_SYSROOT}>>:-Xcompiler=-isysroot;-Xcompiler=${CMAKE_OSX_SYSROOT}>
    # configurations for LLVM
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>,$<CONFIG:Debug>>:-g>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>,$<CONFIG:Release>>:-O3 -DNDEBUG>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>,$<CONFIG:MinSizeRel>>:-O2 -DNDEBUG>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>,$<CONFIG:RelWithDebInfo>>:-O2 -g -DNDEBUG>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>,$<CONFIG:Coverage>>:-Og -g>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>,$<CONFIG:RelWithAssert>>:-O3 -g>
)
# cmake-format: on

#
# Code coverage
#

if(ESPRESSO_BUILD_WITH_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
    # cmake-format: off
    target_compile_options(
      espresso_compiler_flags
      INTERFACE
        $<$<COMPILE_LANGUAGE:CXX,CUDA>:-fprofile-instr-generate -fcoverage-mapping>
    )
    # cmake-format: on
  else()
    # cmake-format: off
    target_compile_options(
      espresso_compiler_flags
      INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:--coverage -fprofile-abs-path>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=--coverage,-fprofile-abs-path>
        # workaround for https://github.com/espressomd/espresso/issues/4943
        $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<BOOL:${ESPRESSO_BUILD_WITH_CCACHE}>>:--coverage -fprofile-abs-path>
    )
    # cmake-format: on
    target_link_libraries(espresso_compiler_flags INTERFACE gcov)
  endif()
endif()

#
# Portability options
#

# prevent 80-bit arithmetic on old Intel processors
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SIZEOF_VOID_P EQUAL 4
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "[xX]86")
  target_compile_options(espresso_compiler_flags
                         INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-ffloat-store>)
endif()

# allow infinities
target_compile_options(
  espresso_compiler_flags
  INTERFACE
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-fno-finite-math-only>
)

#
# Sanitizers
#

if(ESPRESSO_BUILD_WITH_ASAN AND ESPRESSO_BUILD_WITH_MSAN)
  message(
    FATAL_ERROR
      "Address sanitizer and memory sanitizer cannot be enabled simultaneously")
endif()
if(ESPRESSO_BUILD_WITH_ASAN)
  target_compile_options(
    espresso_compiler_flags
    INTERFACE
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:${ESPRESSO_XCOMPILER}-fsanitize=address>
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:${ESPRESSO_XCOMPILER}-g>
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:$<$<NOT:$<CONFIG:Debug>>:${ESPRESSO_XCOMPILER}-O1>>
      $<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer>)
  target_link_libraries(
    espresso_compiler_flags INTERFACE -fsanitize=address -g
                                      $<$<NOT:$<CONFIG:Debug>>:-O1>)
endif()
if(ESPRESSO_BUILD_WITH_MSAN)
  target_compile_options(
    espresso_compiler_flags
    INTERFACE
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:${ESPRESSO_XCOMPILER}-fsanitize=memory>
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:${ESPRESSO_XCOMPILER}-g>
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:$<$<NOT:$<CONFIG:Debug>>:${ESPRESSO_XCOMPILER}-O1>>
      $<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer>)
  target_link_libraries(
    espresso_compiler_flags INTERFACE -fsanitize=memory -g
                                      $<$<NOT:$<CONFIG:Debug>>:-O1>)
endif()
if(ESPRESSO_BUILD_WITH_UBSAN)
  target_compile_options(
    espresso_compiler_flags
    INTERFACE
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:-fsanitize=undefined>
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:${ESPRESSO_XCOMPILER}-g>
      $<$<COMPILE_LANGUAGE:CXX,CUDA>:$<$<NOT:$<CONFIG:Debug>>:${ESPRESSO_XCOMPILER}-O1>>
  )
  target_link_libraries(
    espresso_compiler_flags INTERFACE -fsanitize=undefined -g
                                      $<$<NOT:$<CONFIG:Debug>>:-O1>)
endif()

#
# Static analysis
#

if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
  find_package(ClangTidy "${CMAKE_CXX_COMPILER_VERSION}" EXACT REQUIRED)
  set(ESPRESSO_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  set(ESPRESSO_CUDA_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=--cuda-host-only")
  block(SCOPE_FOR VARIABLES)
  if(ESPRESSO_BUILD_WITH_CUDA)
    # silence casts in cuda_runtime.h (for both C++ and CUDA source files)
    list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-casting-through-void")
    # silence nullptr dereference in cuda::thrust
    list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA
         "-clang-analyzer-core.NonNullParamChecker")
  endif()
  foreach(LANG CXX CUDA)
    espresso_override_clang_tidy_checks(
      ESPRESSO_${LANG}_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
      "${SKIP_CLANG_TIDY_CHECKS_${LANG}}")
    set(ESPRESSO_${LANG}_CLANG_TIDY "${ESPRESSO_${LANG}_CLANG_TIDY}"
        PARENT_SCOPE)
  endforeach()
  endblock()
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

set(ESPRESSO_OPENMP_COMPONENTS CXX)
if(ESPRESSO_BUILD_WITH_CUDA AND NOT CMAKE_CUDA_COMPILER_ID STREQUAL NVHPC)
  # FindOpenMP doesn't fully support NVHPC+CUDA as of CMake 4.3, see ticket:
  # https://gitlab.kitware.com/cmake/cmake/-/work_items/23003
  list(APPEND ESPRESSO_OPENMP_COMPONENTS CUDA)
endif()
# ESPResSo and Kokkos require some of the new features from OpenMP 5.0, but a
# version check cannot be enforced here, since GCC versions 12 to 15 report 4.5
# support for the CXX component, and NVHPC versions <= 26.3 report 4.5 support
# for the CUDA component
find_package(OpenMP REQUIRED COMPONENTS ${ESPRESSO_OPENMP_COMPONENTS})

if(ESPRESSO_BUILD_WITH_FFTW)
  find_package(fftw3 REQUIRED COMPONENTS omp)

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
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_MakeAvailable(heffte)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    foreach(HEFFTE_MODULE IN ITEMS FFTW CUDA)
      if(TARGET Heffte::${HEFFTE_MODULE})
        set(Heffte_${HEFFTE_MODULE}_FOUND true)
      endif()
    endforeach()
    add_library(Heffte::Heffte INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Heffte::Heffte INTERFACE Heffte)
    set_property(TARGET Heffte PROPERTY INSTALL_RPATH "")
    install(TARGETS Heffte LIBRARY DESTINATION "${ESPRESSO_INSTALL_LIBDIR}")
  endif()
endif()

if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/kokkos-src)
  find_package(Kokkos 4.6 QUIET)
endif()
if(NOT DEFINED Kokkos_FOUND OR NOT ${Kokkos_FOUND})
  # cmake-format: off
  FetchContent_Declare(
    kokkos
    GIT_REPOSITORY https://github.com/kokkos/kokkos.git
    GIT_TAG        5.0.2
    OVERRIDE_FIND_PACKAGE
  )
  # cmake-format: on
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
  set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
  set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OFF CACHE BOOL "")
  set(Kokkos_ENABLE_IMPL_VIEW_LEGACY ON CACHE BOOL "")
  set(Kokkos_ENABLE_COMPLEX_ALIGN ON CACHE BOOL "")
  set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON CACHE BOOL "")
  set(Kokkos_ENABLE_HWLOC ON CACHE BOOL "")
  set(Kokkos_ENABLE_DEPRECATED_CODE_4 OFF CACHE BOOL "")
  set(Kokkos_ENABLE_DEPRECATED_CODE_5 OFF CACHE BOOL "")
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
  target_compile_options(
    kokkoscore
    PRIVATE
      $<$<COMPILE_LANG_AND_ID:CXX,NVHPC>:--diag_suppress=initialization_not_reachable>
  )
  # mark kokkos headers as system headers to disable compiler diagnostics
  set_property(
    TARGET kokkos APPEND
    PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
             $<TARGET_PROPERTY:kokkos,INTERFACE_INCLUDE_DIRECTORIES>)
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
  # mark Cabana headers as system headers to disable compiler diagnostics
  set_property(
    TARGET Core APPEND
    PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
             $<TARGET_PROPERTY:Core,INTERFACE_INCLUDE_DIRECTORIES>)
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
  if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/HighFive-src)
    find_package(HighFive 3.3.0 QUIET)
  endif()
  if(NOT DEFINED HighFive_FOUND OR NOT ${HighFive_FOUND})
    # cmake-format: off
    FetchContent_Declare(
      HighFive
      GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
      GIT_TAG        v3.3.0
    )
    # cmake-format: on
    FetchContent_MakeAvailable(HighFive)
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

if(ESPRESSO_BUILD_WITH_NLOPT)
  if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/nlopt-src)
    find_package(NLopt 2.7.1 QUIET)
  endif()
  if(NOT DEFINED NLopt_FOUND OR NOT ${NLopt_FOUND})
    # cmake-format: off
    FetchContent_Declare(
      nlopt
      GIT_REPOSITORY https://github.com/stevengj/nlopt.git
      GIT_TAG        v2.11.0
    )
    # cmake-format: on
    set(BUILD_SHARED_LIBS off)
    set(NLOPT_FORTRAN off CACHE BOOL "")
    set(NLOPT_PYTHON off CACHE BOOL "")
    set(NLOPT_OCTAVE off CACHE BOOL "")
    set(NLOPT_MATLAB off CACHE BOOL "")
    set(NLOPT_GUILE off CACHE BOOL "")
    set(NLOPT_JAVA off CACHE BOOL "")
    set(NLOPT_SWIG off CACHE BOOL "")
    set(NLOPT_LUKSAN off CACHE BOOL "")
    FetchContent_MakeAvailable(nlopt)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    target_compile_options(
      nlopt
      PRIVATE
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=code_is_unreachable>
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=integer_sign_change>
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=deprecated_entity>
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=mixed_enum_type>)
  endif()
endif()

if(ESPRESSO_BUILD_WITH_VALGRIND)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(VALGRIND valgrind REQUIRED)
  if(VALGRIND_FOUND)
    message(STATUS "Found valgrind: ${VALGRIND_INCLUDE_DIRS}")
  endif()
endif()

if(Kokkos_ENABLE_CUDA AND NOT CMAKE_CXX_COMPILER_ID STREQUAL Clang)
  # allow constexpr __host__ function calls in __device__ functions
  target_compile_options(espresso_kokkos_compiler_flags
                         INTERFACE -expt-relaxed-constexpr)
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

# MPIEXEC_PREFLAGS overrides
set(ESPRESSO_MPIEXEC_PREFLAGS "")
# Open MPI 4.x has a bug on NUMA archs that prevents running in singleton mode
set(ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA OFF)
set(ESPRESSO_CPU_MODEL_NAME_OMPI_SINGLETON_NUMA_PATTERN "AMD (EPYC|Ryzen)")

if("${ESPRESSO_MPIEXEC_VENDOR}" STREQUAL "OpenMPI")
  # OpenMPI checks the number of processes against the number of physical cores
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
  set(ESPRESSO_TEST_NT "4"
      CACHE STRING "Maximal number of OpenMP threads to use per test")
  add_library(espresso_tests_compiler_flags INTERFACE)
  add_library(espresso::tests::compiler_flags ALIAS
              espresso_tests_compiler_flags)
  target_compile_options(
    espresso_tests_compiler_flags
    INTERFACE $<$<NOT:$<CXX_COMPILER_ID:NVHPC>>:-Wno-unused-macros>)
  if(ESPRESSO_BUILD_WITH_COVERAGE)
    target_compile_options(
      espresso_tests_compiler_flags
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
  if(NOT EXISTS ${FETCHCONTENT_BASE_DIR}/walberla-src)
    find_package(waLBerla 7.2 QUIET)
    if(waLBerla_FOUND AND NOT TARGET walberla::core AND DEFINED
                                                        walberla_SOURCE_DIR)
      add_subdirectory(${walberla_SOURCE_DIR} ${walberla_BINARY_DIR}
                       EXCLUDE_FROM_ALL)
    endif()
  endif()
  if(NOT DEFINED waLBerla_FOUND OR NOT ${waLBerla_FOUND})
    # cmake-format: off
    FetchContent_Declare(
      walberla
      GIT_REPOSITORY https://i10git.cs.fau.de/walberla/walberla.git
      GIT_TAG        3247aa73 # v7.2 with patches
    )
    # cmake-format: on
    string(REGEX REPLACE "([/\\]walberla)-src$" "\\1-build" walberla_BINARY_DIR
                         "${walberla_SOURCE_DIR}")
    set(WALBERLA_BUILD_TESTS off CACHE BOOL "")
    set(WALBERLA_BUILD_TOOLS off CACHE BOOL "")
    set(WALBERLA_BUILD_LEGACY_LBM off CACHE BOOL "")
    set(WALBERLA_BUILD_BENCHMARKS off CACHE BOOL "")
    set(WALBERLA_BUILD_TUTORIALS off CACHE BOOL "")
    set(WALBERLA_BUILD_SHOWCASES off CACHE BOOL "")
    set(WALBERLA_BUILD_EXAMPLES off CACHE BOOL "")
    set(WALBERLA_BUILD_DOC off CACHE BOOL "")
    set(WALBERLA_BUILD_WITH_PYTHON off CACHE BOOL "")
    set(WALBERLA_LOGLEVEL "WARNING" CACHE STRING "")
    if(ESPRESSO_BUILD_WITH_CUDA)
      set(WALBERLA_BUILD_WITH_CUDA "on" CACHE BOOL "")
      if(NOT ESPRESSO_CUDA_COMPILER STREQUAL "clang")
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
          message(FATAL_ERROR "variable CMAKE_CUDA_ARCHITECTURES is undefined")
        endif()
      endif()
    endif()
    set(WALBERLA_BUILD_WITH_FFTW off CACHE BOOL "")
    set(WALBERLA_BUILD_WITH_OPENMP on CACHE BOOL "")
    set(WALBERLA_BUILD_WITH_FASTMATH off CACHE BOOL "")
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_MakeAvailable(walberla)
    set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
    set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
    add_library(espresso_walberla_sqlite_compiler_flags INTERFACE)
    target_compile_options(
      espresso_walberla_sqlite_compiler_flags
      INTERFACE
        $<$<COMPILE_LANG_AND_ID:C,GNU>:-Wno-discarded-qualifiers>
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,CrayClang,AppleClang,IntelLLVM>:-Wno-misleading-indentation>
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Wno-misleading-indentation>
        -Wno-unused-but-set-variable)
    target_link_libraries(walberla_sqlite
                          PRIVATE espresso_walberla_sqlite_compiler_flags)
    target_compile_options(
      sqlite3
      PRIVATE
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=code_is_unreachable>
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=cast_to_qualified_type>
        $<$<COMPILE_LANG_AND_ID:C,NVHPC>:--diag_suppress=set_but_not_used>)
  endif()
  add_library(espresso_walberla_deps INTERFACE)
  add_library(espresso::walberla_deps ALIAS espresso_walberla_deps)
  target_link_libraries(
    espresso_walberla_deps
    INTERFACE walberla::core walberla::domain_decomposition
              walberla::blockforest walberla::communication walberla::field
              walberla::stencil walberla::vtk
              $<$<BOOL:${WALBERLA_BUILD_WITH_FFTW}>:walberla::fft>
              $<$<BOOL:${WALBERLA_BUILD_WITH_CUDA}>:walberla::gpu>)
  if(ESPRESSO_BUILD_WITH_WALBERLA_AVX)
    function(espresso_avx_flags_callback COMPILER_AVX2_FLAG)
      target_compile_options(
        espresso_avx_flags
        INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${COMPILER_AVX2_FLAG}>
                  -DESPRESSO_BUILD_WITH_AVX_KERNELS)
    endfunction()
    espresso_enable_avx2_support(espresso_avx_flags_callback)
  endif()
endif()

if(ESPRESSO_BUILD_WITH_CALIPER)
  # cmake-format: off
  FetchContent_Declare(
    caliper
    GIT_REPOSITORY https://github.com/LLNL/Caliper.git
    GIT_TAG        v2.14.0
  )
  # cmake-format: on
  set(CALIPER_OPTION_PREFIX on CACHE BOOL "")
  set(CALIPER_WITH_MPI on CACHE BOOL "")
  set(CALIPER_WITH_NVTX ${ESPRESSO_BUILD_WITH_CUDA} CACHE BOOL "")
  set(CALIPER_WITH_CUPTI ${ESPRESSO_BUILD_WITH_CUDA} CACHE BOOL "")
  set(CALIPER_INSTALL_CONFIG off CACHE BOOL "")
  set(CALIPER_INSTALL_HEADERS off CACHE BOOL "")
  set(BUILD_SHARED_LIBS ON)
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_INSTALL_LIBDIR}")
  set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
  FetchContent_MakeAvailable(caliper)
  set(BUILD_SHARED_LIBS ${ESPRESSO_BUILD_SHARED_LIBS_DEFAULT})
  set(CMAKE_INSTALL_LIBDIR "${ESPRESSO_OLD_CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_SHARED_LIBRARY_PREFIX "${ESPRESSO_SHARED_LIBRARY_PREFIX}")
  add_library(espresso_caliper_compiler_flags INTERFACE)
  target_compile_options(
    espresso_caliper_compiler_flags
    INTERFACE $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-free-nonheap-object
              -Wno-deprecated-enum-enum-conversion -Wno-volatile
              -Wno-maybe-uninitialized>)
  target_link_libraries(caliper-common PRIVATE espresso_caliper_compiler_flags)
  target_link_libraries(caliper-runtime PRIVATE espresso_caliper_compiler_flags)
endif()

#
# Set source file properties
#

function(espresso_set_common_target_properties)
  set(TARGET_NAME "${ARGV0}")
  get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)
  if(NOT TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
    # set non-transitive properties; transitive properties should be propagated
    # to dependent targets via INTERFACE targets, e.g. espresso::compiler_flags
    set_target_properties(${TARGET_NAME} PROPERTIES C_EXTENSIONS OFF)
    set_target_properties(${TARGET_NAME} PROPERTIES CXX_EXTENSIONS OFF)
    set_target_properties(${TARGET_NAME} PROPERTIES CUDA_EXTENSIONS OFF)
    if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
      foreach(LANG CXX CUDA)
        set_target_properties(
          ${TARGET_NAME} PROPERTIES ${LANG}_CLANG_TIDY
                                    "${ESPRESSO_${LANG}_CLANG_TIDY}")
      endforeach()
    endif()
    if(ESPRESSO_BUILD_WITH_CUDA)
      if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
        set_target_properties(${TARGET_NAME}
                              PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
      elseif(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
        set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE "CXX")
      endif()
      # handle builds with different GCC versions for C++ and CUDA
      # (https://github.com/espressomd/espresso/issues/5256)
      if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA" AND CMAKE_CXX_COMPILER_ID
                                                      STREQUAL "GNU"
         AND NOT "$ENV{CUDAHOSTCXX}" STREQUAL "$ENV{CXX}"
         AND "$ENV{CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES_EXCLUDE}" MATCHES
             "^/usr/lib/gcc/[^/]+/1[4-9]\$")
        target_link_directories(
          ${TARGET_NAME} BEFORE PUBLIC
          "$ENV{CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES_EXCLUDE}")
      endif()
    endif()
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
