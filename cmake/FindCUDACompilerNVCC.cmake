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

set(CMAKE_CUDA_COMPILER ${CUDA_NVCC_EXECUTABLE})
set(CUDA 1)

execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --version
                OUTPUT_VARIABLE NVCC_VERSION_STRING)

string(REGEX
       REPLACE "^.*Cuda compilation tools, release [0-9\.]+, V([0-9\.]+).*\$"
               "\\1"
               CMAKE_CUDA_COMPILER_VERSION
               "${NVCC_VERSION_STRING}")

if(NOT CUDA_NVCC_EXECUTABLE STREQUAL "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc")
  get_filename_component(NVCC_EXECUTABLE_DIRNAME "${CUDA_NVCC_EXECUTABLE}" DIRECTORY)
  get_filename_component(NVCC_EXECUTABLE_DIRNAME "${NVCC_EXECUTABLE_DIRNAME}" DIRECTORY)
  message(WARNING "Your nvcc (${CUDA_NVCC_EXECUTABLE}) does not appear to match your CUDA libraries (in ${CUDA_TOOLKIT_ROOT_DIR}). While ESPResSo will still compile, you might get unexpected crashes. Please point CUDA_TOOLKIT_ROOT_DIR to your CUDA toolkit path, e.g. by adding -DCUDA_TOOLKIT_ROOT_DIR='${NVCC_EXECUTABLE_DIRNAME}' to your cmake command.")
endif()

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG} -g")
set(CUDA_NVCC_FLAGS_RELEASE "${CUDA_NVCC_FLAGS_RELEASE} -O3 -Xptxas=-O3 -Xcompiler=-O3 -DNDEBUG")
set(CUDA_NVCC_FLAGS_MINSIZEREL "${CUDA_NVCC_FLAGS_MINSIZEREL} -O2 -Xptxas=-O2 -Xcompiler=-Os -DNDEBUG")
set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO} -O2 -g -Xptxas=-O2 -Xcompiler=-O2,-g -DNDEBUG")
set(CUDA_NVCC_FLAGS_COVERAGE "${CUDA_NVCC_FLAGS_COVERAGE} -O3 -g -Xptxas=-O3 -Xcompiler=-Og,-g")
set(CUDA_NVCC_FLAGS_RELWITHASSERT "${CUDA_NVCC_FLAGS_RELWITHASSERT} -O3 -g -Xptxas=-O3 -Xcompiler=-O3,-g")

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_52,code=compute_52 -std=c++${CMAKE_CUDA_STANDARD}")
if(WARNINGS_ARE_ERRORS)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler=-Werror -Xptxas=-Werror")
endif()
if (CMAKE_OSX_SYSROOT)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler=-isysroot -Xcompiler=${CMAKE_OSX_SYSROOT}")
endif()

function(add_gpu_library)
  cuda_add_library(${ARGV})
  set_property(TARGET ${ARGV0} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(${ARGV0} PRIVATE ${CUDA_CUFFT_LIBRARIES})
endfunction()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( CudaCompilerNVCC
                                   REQUIRED_VARS CMAKE_CUDA_COMPILER
                                   VERSION_VAR CMAKE_CUDA_COMPILER_VERSION )
