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

# Verify the HPC SDK C++ and CUDA compilers match,
# include the toolkit libraries and declare a custom
# `add_library()` wrapper function named `espresso_add_gpu_library()`.

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_CUDA_COMPILER_ID)
  message(
    FATAL_ERROR
      "To compile CUDA code with ${CMAKE_CUDA_COMPILER_ID}, the C++ compiler "
      "must be ${CMAKE_CUDA_COMPILER_ID}, not ${CMAKE_CXX_COMPILER_ID}."
  )
endif()

block(PROPAGATE ESPRESSO_NVHPC_VERBOSE_OUTPUT)
  separate_arguments(ESPRESSO_CMAKE_CUDA_FLAGS_LIST NATIVE_COMMAND "${CMAKE_CUDA_FLAGS}")
  execute_process(COMMAND ${CMAKE_CUDA_COMPILER} ${ESPRESSO_CMAKE_CUDA_FLAGS_LIST} ${ARGV} -v
                  ERROR_VARIABLE OUTPUT_A)
  execute_process(COMMAND ${CMAKE_CUDA_COMPILER} ${ESPRESSO_CMAKE_CUDA_FLAGS_LIST} ${ARGV} -V
                  OUTPUT_VARIABLE OUTPUT_B)
  set(ESPRESSO_NVHPC_VERBOSE_OUTPUT "${OUTPUT_A}\n${OUTPUT_B}")
endblock()

string(REGEX MATCH "[Ee]xport NVHPC_CURRENT_CUDA_HOME=([^\n]+)" _ "${ESPRESSO_NVHPC_VERBOSE_OUTPUT}")
set(ESPRESSO_NVHPC_DETECTED_CUDA_DIR ${CMAKE_MATCH_1})
string(REGEX MATCH "[Ee]xport NVHPC_CURRENT_CUDA_VERSION=([^\n]+)" _ "${ESPRESSO_NVHPC_VERBOSE_OUTPUT}")
set(ESPRESSO_NVHPC_DETECTED_CUDA_VERSION ${CMAKE_MATCH_1})
if ("${CMAKE_CUDA_COMPILER_VERSION}" STREQUAL "")
  string(REGEX MATCH "nvc\\+\\+ ([0-9.]+)" _ "${ESPRESSO_NVHPC_VERBOSE_OUTPUT}")
  set(CMAKE_CUDA_COMPILER_VERSION ${CMAKE_MATCH_1})
endif()

if(NOT EXISTS ${ESPRESSO_NVHPC_DETECTED_CUDA_DIR})
  message(FATAL_ERROR "${CMAKE_CUDA_COMPILER_ID} could not automatically detect a compatible CUDA toolkit library")
endif()

if("${CUDAToolkit_ROOT}" STREQUAL "")
  set(CUDAToolkit_ROOT ${ESPRESSO_NVHPC_DETECTED_CUDA_DIR})
elseif(NOT "${ESPRESSO_NVHPC_DETECTED_CUDA_DIR}" STREQUAL "${CUDAToolkit_ROOT}")
  message(
    WARNING
      "${CMAKE_CUDA_COMPILER_ID} CUDA toolkit directory (${ESPRESSO_NVHPC_DETECTED_CUDA_DIR}) "
      "and NVIDIA CUDA toolkit directory (${CUDAToolkit_ROOT}) don't match; try hinting it with "
      "'-D CUDAToolkit_ROOT=\"${ESPRESSO_NVHPC_DETECTED_CUDA_DIR}\"'.")
endif()

if("${CUDAToolkit_VERSION}" STREQUAL "")
  set(CUDAToolkit_VERSION ${ESPRESSO_NVHPC_DETECTED_CUDA_VERSION})
elseif(NOT "${ESPRESSO_NVHPC_DETECTED_CUDA_VERSION}" STREQUAL "${CUDAToolkit_VERSION}")
  message(
    WARNING
      "The nvc++ CUDA toolkit version (${ESPRESSO_NVHPC_DETECTED_CUDA_VERSION}) "
      "and the NVIDIA CUDA toolkit version (${CUDAToolkit_VERSION}) do not match")
endif()

message(
  STATUS
    "Found CUDA toolkit installation: ${ESPRESSO_NVHPC_DETECTED_CUDA_DIR} "
    "(identified by ${CMAKE_CUDA_COMPILER_ID} as CUDA ${ESPRESSO_NVHPC_DETECTED_CUDA_VERSION})")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerNVHPC REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
