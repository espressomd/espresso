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

if(GIT_EXECUTABLE)
  # Get the name of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get branch status
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD --
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE GIT_DIFF_INDEX_RESULT
    OUTPUT_VARIABLE GIT_DIFF_INDEX_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(GIT_DIFF_INDEX_RESULT EQUAL 0)
    set(GIT_STATE "CLEAN")
  else()
    set(GIT_STATE "DIRTY")
  endif()

endif(GIT_EXECUTABLE)

configure_file(${PROJECT_SOURCE_DIR}/src/config/include/config/version.hpp.in
        ${CMAKE_BINARY_DIR}/include/config/version.hpp.tmp)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/include/config/version.hpp.tmp
        ${CMAKE_BINARY_DIR}/include/config/version.hpp)
