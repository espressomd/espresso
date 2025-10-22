#
# Copyright (C) 2017-2022 The ESPResSo project
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

# Find Clang-tidy

# get Clang version
string(REGEX
        REPLACE "^([0-9]+)\\.[0-9]+.*$"
        "\\1"
        CLANG_MAJOR_VERSION
        "${CMAKE_CXX_COMPILER_VERSION}")
string(REGEX
       REPLACE "^[0-9]+\\.([0-9]+).*$"
               "\\1"
               CLANG_MINOR_VERSION
               "${CMAKE_CXX_COMPILER_VERSION}")
# find Clang-tidy
find_program(CLANG_TIDY_EXE
             NAMES "clang-tidy-${CLANG_MAJOR_VERSION}.${CLANG_MINOR_VERSION}" "clang-tidy-${CLANG_MAJOR_VERSION}" "clang-tidy"
             DOC "Path to clang-tidy executable")
if(CLANG_TIDY_EXE)
  execute_process(COMMAND ${CLANG_TIDY_EXE} --version
                OUTPUT_VARIABLE CLANG_TIDY_OUTPUT RESULT_VARIABLE CLANG_TIDY_STATUS ERROR_QUIET)
  if(CLANG_TIDY_STATUS EQUAL 0)
    string(REGEX MATCH "LLVM version ([0-9]+\\.[0-9\\.]+)" CLANG_TIDY_VERSION "${CLANG_TIDY_OUTPUT}")
    string(REGEX MATCH "([0-9\\.]+)" CLANG_TIDY_VERSION "${CLANG_TIDY_VERSION}")
  endif()
endif()

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( ClangTidy REQUIRED_VARS CLANG_TIDY_EXE
                                   VERSION_VAR CLANG_TIDY_VERSION)
