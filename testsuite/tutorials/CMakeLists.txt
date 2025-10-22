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

function(TUTORIAL_TEST)
  espresso_importlib_based_test(${ARGV} TYPE tutorial)
endfunction(TUTORIAL_TEST)

configure_file(test_convert.py test_convert.py)
configure_file(${ESPRESSO_IMPORTLIB_WRAPPER}
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py COPYONLY)
espresso_importlib_based_fixture(SRC ${ESPRESSO_IMPORTLIB_WRAPPER_TEST} TYPE
                                 tutorial)
espresso_importlib_based_fixture(SRC test_convert.py TYPE tutorial)

# create a local copy of the tutorials
set(TUTORIALS_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_tutorials")
add_custom_target(
  local_tutorials COMMAND ${CMAKE_COMMAND} -E remove_directory
                          ${TUTORIALS_DIR} # cleanup
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/doc/tutorials
          ${TUTORIALS_DIR} DEPENDS tutorials_python)

tutorial_test(FILE test_error_analysis_part1.py)
tutorial_test(FILE test_error_analysis_part2.py)
tutorial_test(FILE test_lennard_jones.py)
tutorial_test(FILE test_charged_system.py)
tutorial_test(FILE test_langevin_dynamics.py)
tutorial_test(FILE test_polymers.py SUFFIX rouse)
tutorial_test(FILE test_polymers.py SUFFIX zimm)
tutorial_test(FILE test_lattice_boltzmann_poiseuille_flow.py)
tutorial_test(FILE test_lattice_boltzmann_sedimentation.py)
tutorial_test(FILE test_raspberry_electrophoresis.py)
tutorial_test(FILE test_active_matter.py)
tutorial_test(FILE test_electrokinetics.py)
tutorial_test(FILE test_visualization.py)
tutorial_test(FILE test_ferrofluid_1.py)
tutorial_test(FILE test_ferrofluid_2.py)
tutorial_test(FILE test_ferrofluid_3.py)
tutorial_test(FILE test_constant_pH__ideal.py)
tutorial_test(FILE test_electrodes_1.py)
tutorial_test(FILE test_electrodes_2.py)
tutorial_test(FILE test_constant_pH__interactions.py)
tutorial_test(FILE test_widom_insertion.py)
tutorial_test(FILE test_grand_canonical_monte_carlo.py)
tutorial_test(FILE test_boltzmann_inversion.py)

add_custom_target(
  check_tutorials JOB_POOL console
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          ${ESPRESSO_CTEST_ARGS} --output-on-failure)

add_dependencies(check_tutorials pypresso local_tutorials)
