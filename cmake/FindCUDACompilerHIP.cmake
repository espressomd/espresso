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

# Find the HIP compiler, include its libraries and declare a custom
# `add_library()` wrapper function named `add_gpu_library()`.

set(ROCM_HOME "/opt/rocm" CACHE FILEPATH "Path to AMD ROCm")
list(APPEND CMAKE_MODULE_PATH "${ROCM_HOME}/hip/cmake")
find_package(HIP ${CUDACompilerHIP_FIND_VERSION} MODULE REQUIRED)
# patch HCC_PATH environment variable and reload HIP
if(HIP_VERSION VERSION_LESS 3.1)
  set(HCC_PATH "${HIP_ROOT_DIR}")
else()
  set(HCC_PATH "${ROCM_HOME}/hcc")
endif()
find_package(HIP ${CUDACompilerHIP_FIND_VERSION} MODULE REQUIRED)

set(CUDA 1)
set(HIP 1)

list(APPEND HIP_HCC_FLAGS
       -std=c++${CMAKE_CUDA_STANDARD} -pedantic -Wall -Wextra
       -Wno-sign-compare -Wno-unused-function -Wno-unused-variable
       -Wno-unused-parameter -Wno-missing-braces -Wno-gnu-anonymous-struct
       -Wno-nested-anon-types -Wno-gnu-zero-variadic-macro-arguments
       -Wno-c99-designator -Wno-macro-redefined -Wno-duplicate-decl-specifier
       $<$<VERSION_GREATER_EQUAL:${HIP_VERSION},3.3>:-Wno-deprecated-copy>
       $<$<BOOL:${WARNINGS_ARE_ERRORS}>:-Werror>)

list(APPEND HIP_HIPCC_FLAGS_DEBUG -g)
list(APPEND HIP_HIPCC_FLAGS_RELEASE -O3 -DNDEBUG)
list(APPEND HIP_HIPCC_FLAGS_MINSIZEREL -O2 -DNDEBUG)
list(APPEND HIP_HIPCC_FLAGS_RELWITHDEBINFO -O2 -g -DNDEBUG)
list(APPEND HIP_HIPCC_FLAGS_COVERAGE -O3 -g)
list(APPEND HIP_HIPCC_FLAGS_RELWITHASSERT -O3 -g)

find_library(ROCFFT_LIB name "rocfft" PATHS ${ROCM_HOME}/lib)

function(add_gpu_library)
  hip_add_library(${ARGV})
  set(GPU_TARGET_NAME ${ARGV0})
  set_target_properties(${GPU_TARGET_NAME} PROPERTIES LINKER_LANGUAGE HIP)
  target_link_libraries(${GPU_TARGET_NAME} PRIVATE ${ROCFFT_LIB})
  target_include_directories(${GPU_TARGET_NAME} PRIVATE ${HIP_ROOT_DIR}/include ${ROCM_HOME}/include)
endfunction()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CudaCompilerHIP REQUIRED_VARS HIP_HIPCC_EXECUTABLE VERSION_VAR HIP_VERSION)
