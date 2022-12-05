#
# Copyright (C) 2019-2022 The ESPResSo project
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

set(TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER
    ${CMAKE_CURRENT_BINARY_DIR}/test_importlib_wrapper.py)
configure_file(importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)
configure_file(test_importlib_wrapper.py
               ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})

macro(PYTHON_SCRIPTS_TEST)
  cmake_parse_arguments(TEST "" "FILE;SUFFIX;TYPE" "DEPENDENCIES;LABELS"
                        ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_with_${TEST_SUFFIX}")
  endif()
  set(TEST_FILE_CONFIGURED "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py")
  configure_file(${TEST_FILE} ${TEST_FILE_CONFIGURED})
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency})
  endforeach(dependency)
  string(REGEX REPLACE "^test_" "${TEST_TYPE}_" TEST_NAME ${TEST_NAME})
  add_test(${TEST_NAME} ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
           ${TEST_FILE_CONFIGURED})
  set_tests_properties(${TEST_NAME} PROPERTIES FIXTURES_REQUIRED
                                               IMPORTLIB_WRAPPER)
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}")
  if("gpu" IN_LIST TEST_LABELS AND ESPRESSO_BUILD_WITH_CUDA)
    set_tests_properties(${TEST_NAME} PROPERTIES RESOURCE_LOCK GPU)
  endif()
endmacro(PYTHON_SCRIPTS_TEST)

add_subdirectory(benchmarks)
add_subdirectory(samples)
add_subdirectory(tutorials)
add_subdirectory(utils)
