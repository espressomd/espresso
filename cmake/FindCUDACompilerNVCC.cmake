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

# include the toolkit libraries and declare a custom
# `add_library()` wrapper function named `espresso_add_gpu_library()`.

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

target_compile_options(
  espresso_cuda_flags
  INTERFACE
  $<$<CONFIG:Debug>:-g -G>
  $<$<CONFIG:Release>:-Xptxas=-O3 -Xcompiler=-O3 -DNDEBUG>
  $<$<CONFIG:MinSizeRel>:-Xptxas=-O2 -Xcompiler=-Os -DNDEBUG>
  $<$<CONFIG:RelWithDebInfo>:-Xptxas=-O2 -Xcompiler=-O2,-g -DNDEBUG>
  $<$<CONFIG:Coverage>:-Xptxas=-O3 -Xcompiler=-Og,-g,--coverage,-fprofile-abs-path>
  $<$<CONFIG:RelWithAssert>:-Xptxas=-O3 -Xcompiler=-O3,-g>
  $<$<BOOL:${CMAKE_OSX_SYSROOT}>:-Xcompiler=-isysroot;-Xcompiler=${CMAKE_OSX_SYSROOT}>
  # workaround for https://github.com/espressomd/espresso/issues/4943
  $<$<BOOL:${ESPRESSO_BUILD_WITH_CCACHE}>:$<$<CONFIG:Coverage>:--coverage -fprofile-abs-path>>
)

function(espresso_configure_gpu_target)
  set(TARGET_NAME ${ARGV0})
  set_target_properties(${TARGET_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(${TARGET_NAME} PRIVATE espresso::cuda_flags $<$<CONFIG:Coverage>:gcov>)
  espresso_add_cuda_rpaths(${TARGET_NAME})
endfunction()

function(espresso_add_gpu_library)
  add_library(${ARGV})
  espresso_configure_gpu_target(${ARGV0})
endfunction()

function(espresso_add_gpu_executable)
  add_executable(${ARGV})
  espresso_configure_gpu_target(${ARGV0})
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerNVCC REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
