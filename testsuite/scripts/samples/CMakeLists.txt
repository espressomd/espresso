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

function(SAMPLE_TEST)
  python_scripts_test(${ARGV} TYPE sample) # forward arguments
  set(samples_tests ${samples_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(SAMPLE_TEST)

# configure Python module importlib_wrapper.py
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
         ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP
                                                  IMPORTLIB_WRAPPER)
set(samples_tests ${samples_tests} ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER}
    PARENT_SCOPE)
configure_file(../importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)

# create a local copy of the samples
set(SAMPLES_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_samples")
add_custom_target(
  local_samples COMMAND ${CMAKE_COMMAND} -E remove_directory
                        ${SAMPLES_DIR} # cleanup
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/samples
          ${SAMPLES_DIR})

sample_test(FILE test_billiard.py)
sample_test(FILE test_chamber_game.py)
sample_test(FILE test_constraints.py)
sample_test(FILE test_dancing.py SUFFIX cpu_ft)
sample_test(FILE test_dancing.py SUFFIX cpu_fts)
sample_test(FILE test_diffusion_coefficient.py)
sample_test(FILE test_dpd.py)
sample_test(FILE test_drude_bmimpf6.py SUFFIX cpu)
sample_test(FILE test_drude_bmimpf6.py SUFFIX gpu LABELS "gpu")
# sample_test(FILE test_ekboundaries.py LABELS "gpu") # TODO WALBERLA
sample_test(FILE test_electrophoresis.py)
sample_test(FILE test_espresso_logo.py)
sample_test(FILE test_gibbs_ensemble.py)
sample_test(FILE test_grand_canonical.py)
if(HDF5_FOUND)
  sample_test(FILE test_h5md.py)
  sample_test(FILE test_h5md_trajectory.py)
endif()
sample_test(FILE test_lbf.py)
sample_test(FILE test_lb_profile.py)
# sample_test(FILE test_lb_planar_couette.py) # TODO WALBERLA
sample_test(FILE test_lb_circular_couette.py)
sample_test(FILE test_lb_four_roller_mill.py)
sample_test(FILE test_lj_liquid_distribution.py)
sample_test(FILE test_lj_liquid.py)
sample_test(FILE test_lj_liquid_structurefactor.py)
sample_test(FILE test_load_checkpoint.py)
set_tests_properties(sample_load_checkpoint PROPERTIES FIXTURES_REQUIRED
                                                       sample_save_checkpoint)
sample_test(FILE test_minimal-charged-particles.py)
sample_test(FILE test_minimal-diamond.py)
sample_test(FILE test_minimal-polymer.py)
sample_test(FILE test_observables_correlators.py)
sample_test(FILE test_p3m.py SUFFIX cpu)
sample_test(FILE test_p3m.py SUFFIX gpu LABELS "gpu")
sample_test(FILE test_reaction_methods.py SUFFIX constant_pH_ensemble)
sample_test(FILE test_reaction_methods.py SUFFIX reaction_ensemble)
sample_test(FILE test_reaction_ensemble_complex_reaction.py)
sample_test(FILE test_rigid_body.py)
sample_test(FILE test_save_checkpoint.py)
set_tests_properties(sample_save_checkpoint PROPERTIES FIXTURES_SETUP
                                                       sample_save_checkpoint)
sample_test(FILE test_slice_input.py)
sample_test(FILE test_visualization_bonded.py)
sample_test(FILE test_visualization_cellsystem.py)
sample_test(FILE test_visualization_charged.py)
foreach(
  var_constraint_shape
  wall;sphere;ellipsoid;cylinder;spherocylinder;simplepore;slitpore;torus;hollowconicalfrustum
)
  sample_test(FILE test_visualization_constraints.py SUFFIX
              ${var_constraint_shape})
endforeach(var_constraint_shape)
sample_test(FILE test_visualization_interactive.py)
sample_test(FILE test_visualization_lbboundaries.py)
sample_test(FILE test_visualization_ljliquid.py)
sample_test(FILE test_visualization_elc.py)
sample_test(FILE test_visualization_npt.py)
sample_test(FILE test_visualization_poiseuille.py)
sample_test(FILE test_widom_insertion.py)
# sample_test(FILE test_object_in_fluid__motivation.py) # TODO WALBERLA
sample_test(FILE test_immersed_boundary.py)

add_custom_target(
  check_samples
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          ${ESPRESSO_CTEST_ARGS} --output-on-failure)

add_dependencies(check_samples pypresso local_samples)
