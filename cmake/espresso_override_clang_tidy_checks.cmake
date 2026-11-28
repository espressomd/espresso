#
# Copyright (C) 2024-2026 The ESPResSo project
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

# Override Clang-Tidy checks.
#
# This function appends an extra flag "--checks=..." in the parent scope
# variable whose name is passed as first argument. The second and third
# arguments are the general and the language-specific overrides, respectively.
# This way you can append a common set of overrides plus language-specific
# overrides in two variables named MYPROJECT_CXX_CLANG_TIDY and
# MYPROJECT_CUDA_CLANG_TIDY, which are used to set the CXX_CLANG_TIDY and
# CUDA_CLANG_TIDY properties of CMake targets.
#
# Example:
# ```cmake
# include(espresso_override_clang_tidy_checks)
# if(MYPROJECT_BUILD_WITH_CLANG_TIDY)
#   find_package(ClangTidy "${CMAKE_CXX_COMPILER_VERSION}" EXACT REQUIRED)
#   set(MYPROJECT_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
#   set(MYPROJECT_CUDA_CLANG_TIDY "${CLANG_TIDY_EXE};--extra-arg=--cuda-host-only")
#   set(SKIP_CLANG_TIDY_CHECKS "")
#   set(SKIP_CLANG_TIDY_CHECKS_CXX "")
#   set(SKIP_CLANG_TIDY_CHECKS_CUDA "")
#   # silence false positives in code enclosed in `extern "C" { /* ... */ }`
#   list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
#   if (MYPROJECT_BUILD_WITH_CUDA)
#     # silence casts in cuda_runtime.h (for both C++ and CUDA source files)
#     list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-casting-through-void")
#     # silence nullptr dereference in cuda::thrust (only for CUDA files)
#     list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA "-clang-analyzer-core.NonNullParamChecker")
#   endif()
#   espresso_override_clang_tidy_checks(MYPROJECT_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}" "${SKIP_CLANG_TIDY_CHECKS_CXX}")
#   espresso_override_clang_tidy_checks(MYPROJECT_CUDA_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}" "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
#   set_target_properties(myproject_core PROPERTIES CXX_CLANG_TIDY "${MYPROJECT_CXX_CLANG_TIDY}")
#   set_target_properties(myproject_cuda PROPERTIES CUDA_CLANG_TIDY "${MYPROJECT_CUDA_CLANG_TIDY}")
# endif()
# ```

function(espresso_override_clang_tidy_checks)
  set(VARNAME "${ARGV0}")
  set(ORIGINAL_OPTIONS "${${VARNAME}}")
  unset(CHECKS)
  string(REGEX MATCH ";--?checks=([^;]*)" HAS_CHECKS "${ORIGINAL_OPTIONS}")
  if(HAS_CHECKS)
    string(REPLACE "," ";" CHECKS "${CMAKE_MATCH_1}")
    string(REPLACE "${CMAKE_MATCH_0}" "" ORIGINAL_OPTIONS "${ORIGINAL_OPTIONS}")
  endif()
  list(APPEND CHECKS ${ARGV1})
  list(APPEND CHECKS ${ARGV2})
  if(CHECKS)
    list(JOIN CHECKS "," CHECKS_STRING)
    set(${VARNAME} "${ORIGINAL_OPTIONS};--checks=${CHECKS_STRING}" PARENT_SCOPE)
  endif()
endfunction()
