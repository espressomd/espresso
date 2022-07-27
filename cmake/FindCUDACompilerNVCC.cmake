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

# Find the NVCC compiler, include its libraries and declare a custom
# `add_library()` wrapper function named `add_gpu_library()`.

set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})
set(CUDA 1)

execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --version
                OUTPUT_VARIABLE NVCC_VERSION_STRING)

string(REGEX
       REPLACE "^.*Cuda compilation tools, release [0-9\.]+, V([0-9\.]+).*\$"
               "\\1" CMAKE_CUDA_COMPILER_VERSION "${NVCC_VERSION_STRING}")

get_filename_component(CMAKE_CUDA_COMPILER_TOOLKIT "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" REALPATH)
get_filename_component(CMAKE_CUDA_COMPILER_RESOLVED "${CMAKE_CUDA_COMPILER}" REALPATH)
if(NOT "${CMAKE_CUDA_COMPILER_TOOLKIT}" STREQUAL "${CMAKE_CUDA_COMPILER_RESOLVED}"
   AND NOT INSIDE_DOCKER)
  get_filename_component(NVCC_EXECUTABLE_DIRNAME "${CMAKE_CUDA_COMPILER}" DIRECTORY)
  get_filename_component(NVCC_EXECUTABLE_DIRNAME "${NVCC_EXECUTABLE_DIRNAME}" DIRECTORY)
  message(
    WARNING
      "Your nvcc (${CMAKE_CUDA_COMPILER}) does not appear to match your CUDA libraries (in ${CUDA_TOOLKIT_ROOT_DIR}). While ESPResSo will still compile, you might get unexpected crashes. Please point CUDA_TOOLKIT_ROOT_DIR to your CUDA toolkit path, e.g. by adding -DCUDA_TOOLKIT_ROOT_DIR='${NVCC_EXECUTABLE_DIRNAME}' to your cmake command."
  )
endif()

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

set(CUDA_PROPAGATE_HOST_FLAGS OFF)

add_library(Espresso_cuda_flags INTERFACE)
add_library(Espresso::cuda_flags ALIAS Espresso_cuda_flags)
target_compile_options(
  Espresso_cuda_flags
  INTERFACE
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:>
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},Release>:-Xptxas=-O3 -Xcompiler=-O3 -DNDEBUG>
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},MinSizeRel>:-Xptxas=-O2 -Xcompiler=-Os -DNDEBUG>
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},RelWithDebInfo>:-Xptxas=-O2 -Xcompiler=-O2,-g -DNDEBUG>
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},Coverage>:-Xptxas=-O3 -Xcompiler=-Og,-g>
  $<$<STREQUAL:${CMAKE_BUILD_TYPE},RelWithAssert>:-Xptxas=-O3 -Xcompiler=-O3,-g>
  "--compiler-bindir=${CMAKE_CXX_COMPILER}"
  $<$<BOOL:${WARNINGS_ARE_ERRORS}>:-Xcompiler=-Werror;-Xptxas=-Werror>
  $<$<BOOL:${CMAKE_OSX_SYSROOT}>:-Xcompiler=-isysroot;-Xcompiler=${CMAKE_OSX_SYSROOT}>
)

function(add_gpu_library)
  add_library(${ARGV})
  set(GPU_TARGET_NAME ${ARGV0})
  set_property(TARGET ${GPU_TARGET_NAME} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(${GPU_TARGET_NAME} PRIVATE Espresso::cuda_flags)
  list(APPEND cuda_archs 52)
  if(CMAKE_CUDA_COMPILER_VERSION LESS 11)
    list(APPEND cuda_archs 30)
  endif()
  set_target_properties(${GPU_TARGET_NAME} PROPERTIES CUDA_ARCHITECTURES "${cuda_archs}")
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerNVCC REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
