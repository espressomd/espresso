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

# Verify the NVCC compiler matches the NVIDIA toolkit,
# include the toolkit libraries and declare a custom
# `add_library()` wrapper function named `espresso_add_gpu_library()`.

get_filename_component(ESPRESO_CUDAToolkit_ROOT_RESOLVED "${CUDAToolkit_ROOT}/bin/nvcc" REALPATH)
get_filename_component(ESPRESO_CMAKE_CUDA_COMPILER_RESOLVED "${CMAKE_CUDA_COMPILER}" REALPATH)
if(NOT "${ESPRESO_CUDAToolkit_ROOT_RESOLVED}" STREQUAL "${ESPRESO_CMAKE_CUDA_COMPILER_RESOLVED}"
   AND NOT ESPRESSO_INSIDE_DOCKER)
  get_filename_component(ESPRESSO_NVCC_EXECUTABLE_DIRNAME "${CMAKE_CUDA_COMPILER}" DIRECTORY)
  get_filename_component(ESPRESSO_NVCC_EXECUTABLE_DIRNAME "${ESPRESSO_NVCC_EXECUTABLE_DIRNAME}" DIRECTORY)
  message(
    WARNING
      "Your nvcc compiler (${CMAKE_CUDA_COMPILER}) does not appear to match your CUDA toolkit installation (${CUDAToolkit_ROOT}). While ESPResSo will still compile, you might get unexpected crashes. Try hinting it with '-D CUDAToolkit_ROOT=\"${ESPRESSO_NVCC_EXECUTABLE_DIRNAME}\"'."
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
  $<$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>:-Xcompiler=-Werror;-Xptxas=-Werror>
  $<$<BOOL:${CMAKE_OSX_SYSROOT}>:-Xcompiler=-isysroot;-Xcompiler=${CMAKE_OSX_SYSROOT}>
)

function(espresso_add_gpu_library)
  add_library(${ARGV})
  set(TARGET_NAME ${ARGV0})
  set_target_properties(${TARGET_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(${TARGET_NAME} PRIVATE espresso::cuda_flags)
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerNVCC REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
