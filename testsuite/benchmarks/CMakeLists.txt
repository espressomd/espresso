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
  espresso_importlib_based_test(${ARGV} TYPE benchmark)
endfunction(BENCHMARK_TEST)

configure_file(${ESPRESSO_IMPORTLIB_WRAPPER}
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py COPYONLY)
espresso_importlib_based_fixture(SRC ${ESPRESSO_IMPORTLIB_WRAPPER_TEST} TYPE
                                 benchmark)

# create a local copy of the benchmarks
set(BENCHMARKS_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_benchmarks")
add_custom_target(
  local_benchmarks COMMAND ${CMAKE_COMMAND} -E remove_directory
                           ${BENCHMARKS_DIR} # cleanup
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_SOURCE_DIR}/maintainer/benchmarks ${BENCHMARKS_DIR})

benchmark_test(FILE test_lj.py)
benchmark_test(FILE test_lb.py SUFFIX cpu)
benchmark_test(FILE test_lb.py SUFFIX gpu LABELS "gpu")
benchmark_test(FILE test_p3m.py SUFFIX cpu)
benchmark_test(FILE test_p3m.py SUFFIX gpu LABELS "gpu")
benchmark_test(FILE test_ferrofluid.py)
benchmark_test(FILE test_mc_acid_base_reservoir.py)

add_custom_target(
  check_benchmarks
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          ${ESPRESSO_CTEST_ARGS} --output-on-failure)

add_dependencies(check_benchmarks pypresso local_benchmarks)
