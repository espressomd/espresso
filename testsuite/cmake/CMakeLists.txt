#
# Copyright (C) 2018-2022 The ESPResSo project
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

function(CMAKE_TEST)
  cmake_parse_arguments(TEST "" "FILE" "DEPENDENCIES" ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  configure_file(${TEST_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_FILE} @ONLY)
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency}
                   @ONLY)
  endforeach(dependency)
  set(TEST_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_FILE}")
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_FILE})
  set_tests_properties(${TEST_NAME} PROPERTIES RUN_SERIAL ON)
  set(cmake_tests ${cmake_tests} ${TEST_FILE} PARENT_SCOPE)
endfunction(CMAKE_TEST)

set(ESPRESSO_PYTHON_DIR ${CMAKE_INSTALL_PREFIX}/${Python_STDARCH})

cmake_test(FILE test_install.sh DEPENDENCIES BashUnitTests.sh)

add_custom_target(setup_install COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_custom_target(
  check_cmake_install
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT} -C serial
          --output-on-failure)
add_dependencies(check_cmake_install setup_install)
