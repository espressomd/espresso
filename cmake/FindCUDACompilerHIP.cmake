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

set(CMAKE_CUDA_COMPILER ${HIP_HIPCC_EXECUTABLE})
set(CUDA 1)
set(HIP 1)

execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --version
                OUTPUT_VARIABLE HIPCC_VERSION_STRING)

string(REGEX
       REPLACE "^.*HCC [Cc]lang version ([0-9\.]+).*\$"
               "\\1"
               CMAKE_CUDA_COMPILER_VERSION
               "${HIPCC_VERSION_STRING}")

list(APPEND HIP_HCC_FLAGS "-I${HIP_ROOT_DIR}/include -I${ROCM_HOME}/include -Wno-c99-designator -Wno-macro-redefined -Wno-duplicate-decl-specifier -std=c++${CMAKE_CXX_STANDARD}")
list(APPEND HIP_HCC_FLAGS "-pedantic -Wall -Wextra -Wno-sign-compare -Wno-unused-function -Wno-unused-variable -Wno-unused-parameter -Wno-missing-braces -Wno-gnu-anonymous-struct -Wno-nested-anon-types -Wno-gnu-zero-variadic-macro-arguments")
if(WARNINGS_ARE_ERRORS)
  list(APPEND HIP_HCC_FLAGS "-Werror")
endif()

find_library(ROCFFT_LIB name "rocfft" PATHS "${ROCM_HOME}/lib")

function(add_gpu_library)
  hip_add_library(${ARGV})
  set_target_properties(${ARGV0} PROPERTIES LINKER_LANGUAGE HIP)
  target_link_libraries(${ARGV0} PRIVATE "${ROCFFT_LIB}")
endfunction()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( CudaCompilerHIP
                                   REQUIRED_VARS CMAKE_CUDA_COMPILER
                                   VERSION_VAR CMAKE_CUDA_COMPILER_VERSION )
