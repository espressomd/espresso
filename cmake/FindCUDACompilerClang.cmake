# Copyright (C) 2009-2020 The ESPResSo project
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

set(CMAKE_CUDA_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_CUDA_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
set(CUDA 1)

execute_process(COMMAND ${CMAKE_CUDA_COMPILER} ${CMAKE_CXX_FLAGS} --verbose
                ERROR_VARIABLE CUDA_DIR_STRING)
string(REGEX REPLACE "^.*Found CUDA installation: ([^,]+).*\$" "\\1" CUDA_DIR
                     "${CUDA_DIR_STRING}")
string(REGEX REPLACE "^.*Found CUDA installation: .* version ([0-9\.]+|unknown).*\$"
                     "\\1" CUDA_VERSION "${CUDA_DIR_STRING}")

message(STATUS "Found CUDA-capable host compiler: ${CMAKE_CUDA_COMPILER}")
if(NOT CUDA_DIR_STRING MATCHES "Found CUDA installation" OR CUDA_VERSION STREQUAL "unknown")
  message(FATAL_ERROR "Clang found no compatible CUDA library.")
endif()
message(STATUS "Found CUDA version: ${CUDA_VERSION}")
message(STATUS "Found CUDA installation: ${CUDA_DIR}")

if(CUDA_VERSION VERSION_LESS ${MINIMAL_CUDA_VERSION})
  message(
    FATAL_ERROR
      "${CMAKE_CUDA_COMPILER} was built for CUDA ${CUDA_VERSION}: version does not match requirements (CUDA ${MINIMAL_CUDA_VERSION})."
  )
endif()

set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG} -g")
set(CUDA_NVCC_FLAGS_RELEASE "${CUDA_NVCC_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CUDA_NVCC_FLAGS_MINSIZEREL "${CUDA_NVCC_FLAGS_MINSIZEREL} -O2 -DNDEBUG")
set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO} -O2 -g -DNDEBUG")
set(CUDA_NVCC_FLAGS_COVERAGE "${CUDA_NVCC_FLAGS_COVERAGE} -O3 -g")
set(CUDA_NVCC_FLAGS_RELWITHASSERT "${CUDA_NVCC_FLAGS_RELWITHASSERT} -O3 -g")
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
set(gpu_interface_flags "${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} --cuda-gpu-arch=sm_30")
if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.8.9)
  set(gpu_interface_flags "${gpu_interface_flags} --cuda-gpu-arch=sm_52")
endif()

find_library(
  CUDART_LIBRARY NAMES cudart PATHS ${CUDA_DIR}/lib64 ${CUDA_DIR}/lib
                                    /usr/local/nvidia/lib
                                    /usr/lib/x86_64-linux-gnu NO_DEFAULT_PATH)
find_library(
  CUFFT_LIBRARY NAMES cufft PATHS ${CUDA_DIR}/lib64 ${CUDA_DIR}/lib
                                  /usr/local/nvidia/lib
                                  /usr/lib/x86_64-linux-gnu NO_DEFAULT_PATH)

function(add_gpu_library)
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
  target_link_libraries(${TARGET_NAME} PRIVATE ${CUDA_LIBRARY} ${CUDART_LIBRARY} ${CUFFT_LIBRARY})
  foreach(file ${TARGET_SOURCES})
    if(${file} MATCHES "\\.cu$")
      set_source_files_properties (${file} PROPERTY COMPILE_FLAGS ${gpu_interface_flags})
    endif()
  endforeach()
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CudaCompilerClang REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
