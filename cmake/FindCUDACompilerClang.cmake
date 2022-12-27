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

# Verify the Clang compiler matches the NVIDIA toolkit,
# include the toolkit libraries and declare a custom
# `add_library()` wrapper function named `espresso_add_gpu_library()`.

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_CUDA_COMPILER_ID)
  message(
    FATAL_ERROR
      "To compile CUDA code with ${CMAKE_CUDA_COMPILER_ID}, the C++ compiler must be ${CMAKE_CUDA_COMPILER_ID}, not ${CMAKE_CXX_COMPILER_ID}."
  )
endif()

add_library(espresso_cuda_flags INTERFACE)
add_library(espresso::cuda_flags ALIAS espresso_cuda_flags)

function(espresso_detect_clang_cuda_path)
  separate_arguments(ESPRESSO_CMAKE_CUDA_FLAGS_LIST NATIVE_COMMAND "${CMAKE_CUDA_FLAGS}")
  execute_process(COMMAND ${CMAKE_CUDA_COMPILER} ${ESPRESSO_CMAKE_CUDA_FLAGS_LIST} ${ARGV} --verbose
                  ERROR_VARIABLE ESPRESSO_CLANG_VERBOSE_OUTPUT)
  set(ESPRESSO_CLANG_VERBOSE_OUTPUT ${ESPRESSO_CLANG_VERBOSE_OUTPUT} PARENT_SCOPE)
endfunction()

espresso_detect_clang_cuda_path()
if(NOT ESPRESSO_CLANG_VERBOSE_OUTPUT MATCHES "Found CUDA installation")
  unset(ESPRESSO_CLANG_HINT_CUDA_PATHS)
  if(NOT CMAKE_CUDA_FLAGS MATCHES "--cuda-path")
    foreach(ESPRESSO_UNIX_CUDA_PATH ${CUDAToolkit_ROOT} /usr/lib/cuda /usr/local/cuda)
      if(EXISTS ${ESPRESSO_UNIX_CUDA_PATH})
        espresso_detect_clang_cuda_path("--cuda-path=${ESPRESSO_UNIX_CUDA_PATH}")
        if(ESPRESSO_CLANG_VERBOSE_OUTPUT MATCHES "Found CUDA installation")
          list(APPEND ESPRESSO_CLANG_HINT_CUDA_PATHS "${ESPRESSO_UNIX_CUDA_PATH}")
        endif()
      endif()
    endforeach()
  endif()
  set(ESPRESSO_CLANG_HINT_CUDA_PATHS_STR "")
  if(DEFINED ESPRESSO_CLANG_HINT_CUDA_PATHS)
    list(JOIN ESPRESSO_CLANG_HINT_CUDA_PATHS " " ESPRESSO_CLANG_HINT_CUDA_PATHS_STR)
    set(ESPRESSO_CLANG_HINT_CUDA_PATHS_STR " (possible paths: ${ESPRESSO_CLANG_HINT_CUDA_PATHS_STR})")
  endif()
  message(FATAL_ERROR "${CMAKE_CUDA_COMPILER_ID} could not automatically detect a compatible CUDA library; try hinting one with both '-D CUDAToolkit_ROOT=\"...\"' and '-D CMAKE_CUDA_FLAGS=\"--cuda-path=...\"'${ESPRESSO_CLANG_HINT_CUDA_PATHS_STR}.")
endif()

string(REGEX REPLACE "^.*Found CUDA installation: ([^,]+).*\$" "\\1"
                     ESPRESSO_CLANG_DETECTED_CUDA_DIR "${ESPRESSO_CLANG_VERBOSE_OUTPUT}")
string(REGEX REPLACE "^.*Found CUDA installation: .* version ([0-9\.]+|unknown).*\$"
                     "\\1" ESPRESSO_CLANG_DETECTED_CUDA_VERSION "${ESPRESSO_CLANG_VERBOSE_OUTPUT}")

if(NOT "${ESPRESSO_CLANG_DETECTED_CUDA_DIR}" STREQUAL "${CUDAToolkit_ROOT}")
  set(ESPRESSO_CUDA_TOOLKIT_MISMATCH_WARN "${CMAKE_CUDA_COMPILER_ID} CUDA toolkit directory (${ESPRESSO_CLANG_DETECTED_CUDA_DIR}) and NVIDIA CUDA toolkit directory (${CUDAToolkit_ROOT}) don't match")
  if("${CUDAToolkit_ROOT}" STREQUAL "")
    message(WARNING "${ESPRESSO_CUDA_TOOLKIT_MISMATCH_WARN}; try hinting it with '-D CUDAToolkit_ROOT=\"${ESPRESSO_CLANG_DETECTED_CUDA_DIR}\"'.")
  else()
    message(WARNING "${ESPRESSO_CUDA_TOOLKIT_MISMATCH_WARN}; try hinting it with '-D CMAKE_CUDA_FLAGS=\"--cuda-path=${CUDAToolkit_ROOT}\"'.")
  endif()
endif()

if(NOT CMAKE_CUDA_FLAGS MATCHES "-Wno-unknown-cuda-version")
  set(ESPRESSO_CUDA_TOOLKIT_UNKNOWN_WARN "use '-D CMAKE_CUDA_FLAGS=\"-Wno-unknown-cuda-version\"' to override this check")
  if(ESPRESSO_CLANG_DETECTED_CUDA_VERSION STREQUAL "unknown")
    message(FATAL_ERROR "${CMAKE_CUDA_COMPILER_ID} could not detect the version of the CUDA toolkit library; ${ESPRESSO_CUDA_TOOLKIT_UNKNOWN_WARN}")
  elseif(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.0.0" AND
         CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "13.0.0" AND
         ESPRESSO_CLANG_DETECTED_CUDA_VERSION VERSION_LESS "12.0")
    message(WARNING "${CMAKE_CUDA_COMPILER_ID} ${CMAKE_CUDA_COMPILER_VERSION} doesn't natively support CUDA ${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}; ${ESPRESSO_CUDA_TOOLKIT_UNKNOWN_WARN}.")
    target_compile_options(espresso_cuda_flags INTERFACE -Wno-unknown-cuda-version)
  endif()
endif()
message(STATUS "Found CUDA toolkit installation: ${ESPRESSO_CLANG_DETECTED_CUDA_DIR} (recognized by ${CMAKE_CUDA_COMPILER_ID} as CUDA ${ESPRESSO_CLANG_DETECTED_CUDA_VERSION})")

target_compile_options(
  espresso_cuda_flags
  INTERFACE
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Release>:-O3 -DNDEBUG>
  $<$<CONFIG:MinSizeRel>:-O2 -DNDEBUG>
  $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
  $<$<CONFIG:Coverage>:-O3 -g>
  $<$<CONFIG:RelWithAssert>:-O3 -g>
)

function(espresso_add_gpu_library)
  set(options STATIC SHARED MODULE EXCLUDE_FROM_ALL)
  set(oneValueArgs)
  set(multiValueArgs)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  list(GET ARG_UNPARSED_ARGUMENTS 0 TARGET_NAME)
  list(REMOVE_AT ARG_UNPARSED_ARGUMENTS 0)
  set(TARGET_SOURCES ${ARG_UNPARSED_ARGUMENTS})
  set_source_files_properties(${TARGET_SOURCES} PROPERTIES LANGUAGE "CUDA")
  add_library(${ARGV})
  set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE "CXX")
  target_link_libraries(${TARGET_NAME} PRIVATE espresso::cuda_flags)
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerClang REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
