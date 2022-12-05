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

# Find the NVCC compiler, include its libraries and declare a custom
# `add_library()` wrapper function named `add_gpu_library()`.

set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})

execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --version
                OUTPUT_VARIABLE NVCC_VERSION_STRING)

string(REGEX
       REPLACE "^.*Cuda compilation tools, release [0-9\.]+, V([0-9\.]+).*\$"
               "\\1" CMAKE_CUDA_COMPILER_VERSION "${NVCC_VERSION_STRING}")

get_filename_component(CMAKE_CUDA_COMPILER_TOOLKIT "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" REALPATH)
get_filename_component(CMAKE_CUDA_COMPILER_RESOLVED "${CMAKE_CUDA_COMPILER}" REALPATH)
if(NOT "${CMAKE_CUDA_COMPILER_TOOLKIT}" STREQUAL "${CMAKE_CUDA_COMPILER_RESOLVED}"
   AND NOT ESPRESSO_INSIDE_DOCKER)
  get_filename_component(NVCC_EXECUTABLE_DIRNAME "${CMAKE_CUDA_COMPILER}" DIRECTORY)
  get_filename_component(NVCC_EXECUTABLE_DIRNAME "${NVCC_EXECUTABLE_DIRNAME}" DIRECTORY)
  message(
    WARNING
      "Your nvcc (${CMAKE_CUDA_COMPILER}) does not appear to match your CUDA libraries (in ${CUDA_TOOLKIT_ROOT_DIR}). While ESPResSo will still compile, you might get unexpected crashes. Please point CUDA_TOOLKIT_ROOT_DIR to your CUDA toolkit path, e.g. by adding -DCUDA_TOOLKIT_ROOT_DIR='${NVCC_EXECUTABLE_DIRNAME}' to your cmake command."
  )
endif()

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

set(CUDA_PROPAGATE_HOST_FLAGS OFF)

add_library(espresso_cuda_flags INTERFACE)
add_library(espresso::cuda_flags ALIAS espresso_cuda_flags)
target_compile_options(
  espresso_cuda_flags
  INTERFACE
  $<$<CONFIG:Debug>:>
  $<$<CONFIG:Release>:-Xptxas=-O3 -Xcompiler=-O3 -DNDEBUG>
  $<$<CONFIG:MinSizeRel>:-Xptxas=-O2 -Xcompiler=-Os -DNDEBUG>
  $<$<CONFIG:RelWithDebInfo>:-Xptxas=-O2 -Xcompiler=-O2,-g -DNDEBUG>
  $<$<CONFIG:Coverage>:-Xptxas=-O3 -Xcompiler=-Og,-g>
  $<$<CONFIG:RelWithAssert>:-Xptxas=-O3 -Xcompiler=-O3,-g>
  "--compiler-bindir=${CMAKE_CXX_COMPILER}"
  $<$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>:-Xcompiler=-Werror;-Xptxas=-Werror>
  $<$<BOOL:${CMAKE_OSX_SYSROOT}>:-Xcompiler=-isysroot;-Xcompiler=${CMAKE_OSX_SYSROOT}>
)

function(espresso_add_gpu_library)
  add_library(${ARGV})
  set(TARGET_NAME ${ARGV0})
  set_target_properties(${TARGET_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(${TARGET_NAME} PRIVATE espresso::cuda_flags)
  list(APPEND cuda_archs 75) # RTX-2000 series (Turing)
  list(APPEND cuda_archs 61) # GTX-1000 series (Pascal)
  if(CMAKE_CUDA_COMPILER_VERSION LESS 11)
    list(APPEND cuda_archs 52) # GTX-900 series (Maxwell)
  endif()
  set_target_properties(${TARGET_NAME} PROPERTIES CUDA_ARCHITECTURES "${cuda_archs}")
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerNVCC REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
