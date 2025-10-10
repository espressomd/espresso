#
# Copyright (C) 2015-2024 The ESPResSo project
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

if(ESPRESSO_BUILD_WITH_COVERAGE_PYTHON)
  set(PYPRESSO_OPTIONS --coverage)
endif()

set(ESPRESSO_IMPORTLIB_WRAPPER
    ${CMAKE_SOURCE_DIR}/maintainer/parsing/importlib_wrapper.py)
set(ESPRESSO_IMPORTLIB_WRAPPER_TEST
    ${CMAKE_CURRENT_BINARY_DIR}/test_importlib_wrapper.py)
configure_file(test_importlib_wrapper.py test_importlib_wrapper.py)

# setup fixtures for testsuites that dynamically load simulation scripts
macro(ESPRESSO_IMPORTLIB_BASED_FIXTURE)
  cmake_parse_arguments(TEST "" "SRC;TYPE" "" ${ARGN})
  cmake_path(GET TEST_SRC STEM TEST_NAME)
  add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_BINARY_DIR}/pypresso
                                     ${PYPRESSO_OPTIONS} ${TEST_SRC})
  # cmake-format: off
  set_tests_properties(
    ${TEST_NAME} PROPERTIES FIXTURES_SETUP "FIXTURE_${TEST_TYPE}"
                            SKIP_RETURN_CODE 5
                            PROCESSORS 1
                            ENVIRONMENT "OMP_PROC_BIND=false"
                            ENVIRONMENT "OMP_NUM_THREADS=1")
  # cmake-format: on
endmacro(ESPRESSO_IMPORTLIB_BASED_FIXTURE)

# setup tests for testsuites that dynamically load simulation scripts
macro(ESPRESSO_IMPORTLIB_BASED_TEST)
  cmake_parse_arguments(TEST "" "FILE;SUFFIX;TYPE" "DEPENDENCIES;LABELS"
                        ${ARGN})
  if(DEFINED TEST_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Cannot parse arguments ${TEST_UNPARSED_ARGUMENTS}")
  endif()
  cmake_path(GET TEST_FILE STEM TEST_NAME)
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_with_${TEST_SUFFIX}")
  endif()
  set(TEST_FILE_CONFIGURED "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py")
  configure_file(${TEST_FILE} ${TEST_FILE_CONFIGURED})
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency})
  endforeach(dependency)
  string(REGEX REPLACE "^test_" "${TEST_TYPE}_" TEST_NAME ${TEST_NAME})
  add_test(NAME ${TEST_NAME}
           COMMAND ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
                   ${TEST_FILE_CONFIGURED})
  list(APPEND TEST_LABELS "${TYPE}_test")
  # cmake-format: off
  set_tests_properties(
    ${TEST_NAME} PROPERTIES FIXTURES_REQUIRED "FIXTURE_${TEST_TYPE}"
                            LABELS "${TEST_LABELS}"
                            SKIP_RETURN_CODE 5
                            PROCESSORS 1
                            ENVIRONMENT "OMP_PROC_BIND=false"
                            ENVIRONMENT "OMP_NUM_THREADS=1")
  # cmake-format: on
  if("gpu" IN_LIST TEST_LABELS AND ESPRESSO_BUILD_WITH_CUDA)
    set_tests_properties(${TEST_NAME} PROPERTIES RESOURCE_LOCK GPU)
  endif()
endmacro(ESPRESSO_IMPORTLIB_BASED_TEST)

add_subdirectory(python)
add_subdirectory(samples)
add_subdirectory(tutorials)
add_subdirectory(benchmarks)
add_subdirectory(cmake)
add_subdirectory(scripts)
