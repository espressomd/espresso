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
  python_scripts_test(${ARGV} TYPE tutorial) # forward arguments
  set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(TUTORIAL_TEST)

# configure Python modules importlib_wrapper.py and convert.py
configure_file(test_convert.py ${CMAKE_CURRENT_BINARY_DIR})
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
         ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
add_test(convert ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
         ${CMAKE_CURRENT_BINARY_DIR}/test_convert.py)
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP
                                                  IMPORTLIB_WRAPPER)
set_tests_properties(convert PROPERTIES FIXTURES_SETUP IMPORTLIB_WRAPPER)
set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER}
    PARENT_SCOPE)
configure_file(../importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)

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
# TODO WALBERLA: simulation is too slow without GPU
#tutorial_test(FILE test_polymers.py SUFFIX zimm LABELS gpu)
tutorial_test(FILE test_lattice_boltzmann_poiseuille_flow.py)
tutorial_test(FILE test_lattice_boltzmann_sedimentation.py)
# TODO WALBERLA: simulation is too short
# tutorial_test(FILE test_raspberry_electrophoresis.py)
tutorial_test(FILE test_active_matter.py)
# TODO WALBERLA: requires ekin
# tutorial_test(FILE test_electrokinetics.py LABELS gpu)
tutorial_test(FILE test_visualization.py)
tutorial_test(FILE test_ferrofluid_1.py)
tutorial_test(FILE test_ferrofluid_2.py)
tutorial_test(FILE test_ferrofluid_3.py)
tutorial_test(FILE test_constant_pH__ideal.py)
tutorial_test(FILE test_constant_pH__interactions.py)
tutorial_test(FILE test_widom_insertion.py)

add_custom_target(
  check_tutorials
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          ${ESPRESSO_CTEST_ARGS} --output-on-failure)

add_dependencies(check_tutorials pypresso local_tutorials)
