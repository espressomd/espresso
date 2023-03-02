#
# Copyright (C) 2022-2023 The ESPResSo project
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

function(espresso_enable_avx2_support callback)
  set(COMPILER_AVX2_FLAG "")
  foreach(FLAG_NAME "-mavx2" "/arch:AVX2")
    string(REGEX REPLACE "[^0-9A-Za-z_]" "_" FLAG_VARIABLE "${FLAG_NAME}")
    check_cxx_compiler_flag("${flag_name}"
                            COMPILER_HAS_${FLAG_VARIABLE}_FLAG_RESULT)
    if(COMPILER_HAS_${FLAG_VARIABLE}_FLAG_RESULT)
      set(COMPILER_AVX2_FLAG "${FLAG_NAME}")
      cmake_language(CALL ${callback} "${COMPILER_AVX2_FLAG}")
      break()
    endif()
  endforeach()
  if(COMPILER_AVX2_FLAG STREQUAL "")
    message(
      FATAL_ERROR
        "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} doesn't support AVX2-specific compiler flags."
    )
  endif()
  if(NOT COMPILER_AVX2_FLAG STREQUAL "/arch:AVX2")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -march=native -E -v - INPUT_FILE /dev/null
      OUTPUT_VARIABLE MARCH_NATIVE_OUTPUT_STRING
      ERROR_VARIABLE MARCH_NATIVE_OUTPUT_STRING)
    if(NOT "${MARCH_NATIVE_OUTPUT_STRING}" MATCHES "[ \n](\\+avx2|-mavx2|-D__AVX2__)[ \n]")
      message(
        FATAL_ERROR
          "AVX2 not supported on this CPU architecture according to ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. While ESPResSo will still compile, you will trigger SIGILL when calling AVX functions."
      )
    endif()
  endif()
  set(CMAKE_REQUIRED_FLAGS_BACKUP "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${COMPILER_AVX2_FLAG}")
  check_cxx_source_compiles(
    "#include <immintrin.h>
       __m256i xi_i = _mm256_set_epi32(1, 2, 3, 4, 5, 6, 7, 8);
       __m256  xi_s = _mm256_set_ps(0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f);
       __m256d xi_d = _mm256_set_pd(0.0, 1.0, 2.0, 3.0);
       int main() {}
      " COMPILER_HAS_AVX2_SUPPORT)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_BACKUP}")
  if(NOT COMPILER_HAS_AVX2_SUPPORT)
    message(
      FATAL_ERROR
        "Cannot execute a simple AVX2 program compiled by ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}."
    )
  endif()
endfunction()

