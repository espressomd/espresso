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

# Find the Clang compiler, include its libraries and declare a custom
# `add_library()` wrapper function named `add_gpu_library()`.

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(
    FATAL_ERROR
      "To compile CUDA code with Clang, the C++ compiler must be Clang, not ${CMAKE_CXX_COMPILER_ID}."
  )
endif()

add_library(espresso_cuda_flags INTERFACE)
add_library(espresso::cuda_flags ALIAS espresso_cuda_flags)

function(detect_clang_cuda_path)
  execute_process(COMMAND ${CMAKE_CUDA_COMPILER} ${CMAKE_CXX_FLAGS} --verbose
                  ERROR_VARIABLE CLANG_VERBOSE_OUTPUT)
  if(CLANG_VERBOSE_OUTPUT MATCHES "Found CUDA installation")
    set(CLANG_VERBOSE_OUTPUT ${CLANG_VERBOSE_OUTPUT} PARENT_SCOPE)
    return()
  endif()
  if(NOT CMAKE_CXX_FLAGS MATCHES "--cuda-path")
    foreach(unix_cuda_path /usr/lib/cuda /usr/local/cuda)
      if(EXISTS ${unix_cuda_path})
        execute_process(COMMAND ${CMAKE_CUDA_COMPILER} ${CMAKE_CXX_FLAGS}
                        "--cuda-path=${unix_cuda_path}" --verbose
                        ERROR_VARIABLE CLANG_VERBOSE_OUTPUT)
        if(CLANG_VERBOSE_OUTPUT MATCHES "Found CUDA installation")
          set(CLANG_VERBOSE_OUTPUT ${CLANG_VERBOSE_OUTPUT} PARENT_SCOPE)
          message(STATUS "Clang did not automatically detect a compatible CUDA library; adding compiler flag --cuda-path=${unix_cuda_path}")
          target_compile_options(espresso_cuda_flags INTERFACE "--cuda-path=${unix_cuda_path}")
          return()
        endif()
      endif()
    endforeach()
  endif()
endfunction()

set(CMAKE_CUDA_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_CUDA_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})

detect_clang_cuda_path()
string(REGEX REPLACE "^.*Found CUDA installation: ([^,]+).*\$" "\\1" CUDA_DIR
                     "${CLANG_VERBOSE_OUTPUT}")
string(REGEX REPLACE "^.*Found CUDA installation: .* version ([0-9\.]+|unknown).*\$"
                     "\\1" CUDA_VERSION "${CLANG_VERBOSE_OUTPUT}")
message(STATUS "Found CUDA-capable host compiler: ${CMAKE_CUDA_COMPILER}")
if(NOT CLANG_VERBOSE_OUTPUT MATCHES "Found CUDA installation" OR CUDA_VERSION STREQUAL "unknown")
  message(STATUS "Clang did not automatically detect a compatible CUDA library; adding compiler flag -Wno-unknown-cuda-version")
  target_compile_options(espresso_cuda_flags INTERFACE -Wno-unknown-cuda-version)
  message(STATUS "Found CUDA version: ${CUDAToolkit_VERSION}")
  message(STATUS "Found CUDA installation: ${CUDAToolkit_LIBRARY_DIR}")
else()
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.0.0" AND
     CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "13.0.0" AND
     CUDA_VERSION VERSION_GREATER_EQUAL "11.0" AND
     CUDA_VERSION VERSION_LESS "12.0")
    message(STATUS "Clang ${CMAKE_CXX_COMPILER_VERSION} doesn't natively support CUDA ${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}; adding compiler flag -Wno-unknown-cuda-version")
    target_compile_options(espresso_cuda_flags INTERFACE -Wno-unknown-cuda-version)
  endif()
  message(STATUS "Found CUDA version: ${CUDAToolkit_VERSION} (recognized by Clang as ${CUDA_VERSION})")
  message(STATUS "Found CUDA installation: ${CUDA_DIR}")
endif()
set(CUDA_VERSION ${CUDAToolkit_VERSION})

target_compile_options(
  espresso_cuda_flags
  INTERFACE
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Release>:-O3 -DNDEBUG>
  $<$<CONFIG:MinSizeRel>:-O2 -DNDEBUG>
  $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
  $<$<CONFIG:Coverage>:-O3 -g>
  $<$<CONFIG:RelWithAssert>:-O3 -g>
  # GTX-900 series (Maxwell)
  $<$<VERSION_LESS:${CMAKE_CUDA_COMPILER_VERSION},12>:--cuda-gpu-arch=sm_52>
  # GTX-1000 series (Pascal)
  $<$<VERSION_GREATER_EQUAL:${CMAKE_CUDA_COMPILER_VERSION},10>:--cuda-gpu-arch=sm_61>
  # RTX-2000 series (Turing)
  # With Clang 14+, architectures sm_70+ are only supported with Thrust 1.11+
  # from CUDA 11.3+, for details see https://github.com/NVIDIA/cub/pull/170
  $<$<AND:$<VERSION_GREATER_EQUAL:${CMAKE_CUDA_COMPILER_VERSION},10>,$<OR:$<VERSION_LESS:${CMAKE_CUDA_COMPILER_VERSION},14>,$<VERSION_GREATER_EQUAL:${CUDA_VERSION},11.3.0>>>:--cuda-gpu-arch=sm_75>
)

function(espresso_add_gpu_library)
  set(options STATIC SHARED MODULE EXCLUDE_FROM_ALL)
  set(oneValueArgs)
  set(multiValueArgs)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  list(GET ARG_UNPARSED_ARGUMENTS 0 TARGET_NAME)
  list(REMOVE_AT ARG_UNPARSED_ARGUMENTS 0)
  set(TARGET_SOURCES ${ARG_UNPARSED_ARGUMENTS})
  set_source_files_properties(${TARGET_SOURCES} PROPERTIES LANGUAGE "CXX")
  add_library(${ARGV})
  set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE "CXX")
  target_link_libraries(${TARGET_NAME} PRIVATE espresso::cuda_flags)
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerClang REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
