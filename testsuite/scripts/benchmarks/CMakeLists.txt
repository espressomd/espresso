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

function(BENCHMARK_TEST)
  python_scripts_test(${ARGV} TYPE benchmark) # forward arguments
  set(benchmarks_tests ${benchmarks_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(BENCHMARK_TEST)

# configure Python module importlib_wrapper.py
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
         ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP
                                                  IMPORTLIB_WRAPPER)
set(benchmarks_tests ${benchmarks_tests}
                     ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER} PARENT_SCOPE)
configure_file(../importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)

# create a local copy of the benchmarks
set(BENCHMARKS_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_benchmarks")
add_custom_target(
  local_benchmarks COMMAND ${CMAKE_COMMAND} -E remove_directory
                           ${BENCHMARKS_DIR} # cleanup
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_SOURCE_DIR}/maintainer/benchmarks ${BENCHMARKS_DIR})

benchmark_test(FILE test_lj.py)
benchmark_test(FILE test_lb.py)
benchmark_test(FILE test_p3m.py)
benchmark_test(FILE test_ferrofluid.py)
benchmark_test(FILE test_mc_acid_base_reservoir.py)

add_custom_target(
  check_benchmarks
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          ${ESPRESSO_CTEST_ARGS} --output-on-failure)

add_dependencies(check_benchmarks pypresso local_benchmarks)
