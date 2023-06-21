#
# Copyright (C) 2016-2022 The ESPResSo project
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

#[=======================================================================[.rst:

.. command:: python_test

  Configure a python test case.

  Flag arguments:

  * ``NO_MPI``: run test case without MPI, i.e. in singleton mode. Only use
    this option when relevant to the test, e.g. to check signal propagation
    (MPI vendors introduce their own implementation-defined signal handlers).

  Single-value arguments:

  * ``FILE``: testcase filename
  * ``MAX_NUM_PROC``: maximal number of MPI ranks to use. In the CI pipeline,
    there is one CI job that reduces the number of MPI ranks to 3.
  * ``GPU_SLOTS``: number of GPU slots consumed by the test. The total number
    of slots can be found in :file:`resources.json`. One slot corresponds to
    roughly 250 MiB of VRAM. All GPU test cases must specify this value to
    avoid oversubscribing the GPU with the lowest amount of VRAM among all CI
    runners (note that one GPU is often used to run 2 CI jobs in parallel).
  * ``SUFFIX``: suffix to append to the test name. Meant to disambiguate two
    test cases sharing the same file but using different parameters or number
    of MPI ranks.

  Multi-value arguments:

  * ``DEPENDS``: test cases that must pass before this test can run.
  * ``DEPENDENCIES``: file dependencies.
  * ``LABELS``: extra test labels. Labels are automatically generated based
    on the number of GPU slots and MPI ranks. Extra labels include ``long``
    for statistical tests that can safely be skipped in code coverage builds.
  * ``ARGUMENTS``: arguments to pass to the test case command line.

#]=======================================================================]
function(python_test)
  cmake_parse_arguments(TEST "NO_MPI" "FILE;MAX_NUM_PROC;GPU_SLOTS;SUFFIX"
                        "DEPENDS;DEPENDENCIES;LABELS;ARGUMENTS" ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  set(TEST_FILE_CONFIGURED "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py")
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_${TEST_SUFFIX}")
  endif()
  configure_file(${TEST_FILE} ${TEST_FILE_CONFIGURED} COPYONLY)
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency}
                   COPYONLY)
  endforeach(dependency)
  set(TEST_FILE ${TEST_FILE_CONFIGURED})

  if(NOT DEFINED TEST_MAX_NUM_PROC)
    if(${TEST_NO_MPI})
      set(TEST_MAX_NUM_PROC 1)
    else()
      set(TEST_MAX_NUM_PROC 4)
    endif()
  elseif(${TEST_NO_MPI} AND NOT ${TEST_MAX_NUM_PROC} EQUAL 1)
    message(FATAL_ERROR "NO_MPI and MAX_NUM_PROC are mutually exclusive")
  endif()
  if(NOT DEFINED TEST_GPU_SLOTS)
    set(TEST_GPU_SLOTS 0)
  endif()

  if(${TEST_MAX_NUM_PROC} GREATER ${ESPRESSO_TEST_NP})
    set(TEST_NUM_PROC ${ESPRESSO_TEST_NP})
  else()
    set(TEST_NUM_PROC ${TEST_MAX_NUM_PROC})
  endif()

  if(EXISTS ${MPIEXEC} AND NOT ${TEST_NO_MPI})
    espresso_set_mpiexec_tmpdir(${TEST_NAME})
    add_test(
      NAME ${TEST_NAME}
      COMMAND
        ${MPIEXEC} ${ESPRESSO_MPIEXEC_OVERSUBSCRIBE} ${MPIEXEC_NUMPROC_FLAG}
        ${TEST_NUM_PROC} ${MPIEXEC_PREFLAGS} ${ESPRESSO_MPIEXEC_TMPDIR}
        ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
        ${TEST_FILE_CONFIGURED} ${TEST_ARGUMENTS} ${MPIEXEC_POSTFLAGS})
  else()
    add_test(${TEST_NAME} ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
             ${TEST_FILE_CONFIGURED} ${TEST_ARGUMENTS})
  endif()
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${TEST_NUM_PROC}
                                               DEPENDS "${TEST_DEPENDS}")

  if(${TEST_GPU_SLOTS} GREATER 0 AND ESPRESSO_BUILD_WITH_CUDA)
    set_property(TEST ${TEST_NAME} PROPERTY RESOURCE_GROUPS
                                            "gpus:${TEST_GPU_SLOTS}")
  endif()

  if(${TEST_GPU_SLOTS} GREATER 0)
    list(APPEND TEST_LABELS "gpu")
  endif()
  if(${TEST_NUM_PROC} EQUAL 2)
    list(APPEND TEST_LABELS "parallel")
  endif()
  if(${TEST_NUM_PROC} EQUAL 3)
    list(APPEND TEST_LABELS "parallel_odd")
  endif()
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}")

  set(python_tests ${python_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(python_test)

#[=======================================================================[.rst:

.. command:: checkpoint_test

  Configure a python checkpoint test case.

  This is a specialization of python_test() for the checkpointing test cases.
  It uses sane default values that can be overriden if necessary. It has one
  extra single-value argument ``MODES`` that lists the features to activate.
  See the comment lines at the point of use for more details.

#]=======================================================================]
function(checkpoint_test)
  cmake_parse_arguments(TEST "" "MODES;MAX_NUM_PROC;GPU_SLOTS;SUFFIX" "LABELS"
                        ${ARGN})
  if(NOT DEFINED TEST_MAX_NUM_PROC)
    set(TEST_MAX_NUM_PROC 4)
  endif()
  if(NOT DEFINED TEST_GPU_SLOTS)
    set(TEST_GPU_SLOTS 0)
  endif()
  if(TEST_SUFFIX)
    set(TEST_ARGUMENTS "Test_suffix_${TEST_SUFFIX}__${TEST_MODES}")
    set(TEST_SUFFIX "_${TEST_MODES}_${TEST_SUFFIX}")
  else()
    set(TEST_ARGUMENTS "Test__${TEST_MODES}")
    set(TEST_SUFFIX "_${TEST_MODES}")
  endif()
  python_test(
    FILE
    save_checkpoint.py
    MAX_NUM_PROC
    ${TEST_MAX_NUM_PROC}
    LABELS
    ${TEST_LABELS}
    SUFFIX
    ${TEST_SUFFIX}
    ARGUMENTS
    ${TEST_ARGUMENTS}
    GPU_SLOTS
    ${TEST_GPU_SLOTS}
    DEPENDENCIES
    unittest_generator.py)
  python_test(
    FILE
    test_checkpoint.py
    MAX_NUM_PROC
    ${TEST_MAX_NUM_PROC}
    GPU_SLOTS
    ${TEST_GPU_SLOTS}
    LABELS
    ${TEST_LABELS}
    SUFFIX
    ${TEST_SUFFIX}
    ARGUMENTS
    ${TEST_ARGUMENTS}
    DEPENDENCIES
    unittest_generator.py
    DEPENDS
    save_checkpoint_${TEST_SUFFIX})
endfunction(checkpoint_test)

# Checkpoint tests run on 4 cores (can be overriden with MAX_NUM_PROC). The
# combination of modes to activate is stored in MODES. A mode consists of a
# feature with zero or more options; separate features with 2 underscores and
# options with 1 underscore (options can appear in any order). For example,
# "p3m_cpu__lb_cpu_ascii" generates modes P3M, P3M.CPU, LB, LB.CPU, LB.ASCII.
checkpoint_test(MODES therm_lb__p3m_cpu__lj__lb_walberla_cpu_ascii SUFFIX
                1_core MAX_NUM_PROC 1)
checkpoint_test(MODES therm_lb__p3m_cpu__lj__lb_walberla_cpu_ascii)
checkpoint_test(MODES therm_lb__elc_cpu__lj__lb_walberla_cpu_binary)
checkpoint_test(MODES therm_lb__elc_gpu__lj__lb_walberla_cpu_ascii GPU_SLOTS 3)
checkpoint_test(MODES therm_lb__p3m_gpu__lj__lb_walberla_cpu_binary GPU_SLOTS 3)
checkpoint_test(MODES therm_npt__int_npt)
checkpoint_test(MODES int_sd__lj)
checkpoint_test(MODES dp3m_cpu__therm_langevin__int_nvt)
checkpoint_test(MODES therm_dpd__int_nvt)
checkpoint_test(MODES scafacos__therm_bd__int_bd)
checkpoint_test(MODES therm_sdm__int_sdm)

python_test(FILE bond_breakage.py MAX_NUM_PROC 4)
python_test(FILE cell_system.py MAX_NUM_PROC 4)
python_test(FILE get_neighbors.py MAX_NUM_PROC 4)
python_test(FILE get_neighbors.py MAX_NUM_PROC 3 SUFFIX 3_cores)
python_test(FILE tune_skin.py MAX_NUM_PROC 1)
python_test(FILE code_info.py MAX_NUM_PROC 1)
python_test(FILE constraint_homogeneous_magnetic_field.py MAX_NUM_PROC 4)
python_test(FILE cutoffs.py MAX_NUM_PROC 4)
python_test(FILE cutoffs.py MAX_NUM_PROC 1 SUFFIX 1_core)
python_test(FILE constraint_shape_based.py MAX_NUM_PROC 2)
python_test(FILE coulomb_cloud_wall.py MAX_NUM_PROC 4 GPU_SLOTS 3)
python_test(FILE coulomb_tuning.py MAX_NUM_PROC 4 GPU_SLOTS 3 LABELS long)
python_test(FILE accumulator_correlator.py MAX_NUM_PROC 4)
python_test(FILE accumulator_mean_variance.py MAX_NUM_PROC 4)
python_test(FILE accumulator_time_series.py MAX_NUM_PROC 1)
python_test(FILE dawaanr-and-dds-gpu.py MAX_NUM_PROC 1 GPU_SLOTS 1)
python_test(FILE dawaanr-and-bh-gpu.py MAX_NUM_PROC 1 GPU_SLOTS 1)
python_test(FILE dds-and-bh-gpu.py MAX_NUM_PROC 4 GPU_SLOTS 3)
python_test(FILE electrostatic_interactions.py MAX_NUM_PROC 2)
python_test(FILE engine_langevin.py MAX_NUM_PROC 4)
python_test(FILE engine_lb.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE engine_lb.py MAX_NUM_PROC 1 GPU_SLOTS 1 SUFFIX 1_core)
python_test(FILE icc.py MAX_NUM_PROC 4)
python_test(FILE icc_interface.py MAX_NUM_PROC 1 GPU_SLOTS 1)
python_test(FILE mass-and-rinertia_per_particle.py MAX_NUM_PROC 2 LABELS long)
python_test(FILE integrate.py MAX_NUM_PROC 4)
python_test(FILE interactions_bond_angle.py MAX_NUM_PROC 4)
python_test(FILE interactions_bonded_interface.py MAX_NUM_PROC 4)
python_test(FILE interactions_bonded_interface.py MAX_NUM_PROC 1 SUFFIX 1_core)
python_test(FILE interactions_bonded.py MAX_NUM_PROC 2)
python_test(FILE interactions_dihedral.py MAX_NUM_PROC 4)
python_test(FILE interactions_non-bonded_interface.py MAX_NUM_PROC 4)
python_test(FILE interactions_non-bonded.py MAX_NUM_PROC 4)
python_test(FILE observables.py MAX_NUM_PROC 4)
python_test(FILE particle.py MAX_NUM_PROC 4)
python_test(FILE pressure.py MAX_NUM_PROC 4)
python_test(FILE scafacos_dipoles_1d_2d.py MAX_NUM_PROC 4)
python_test(FILE scafacos_interface.py MAX_NUM_PROC 2)
python_test(FILE long_range_actors.py MAX_NUM_PROC 4 GPU_SLOTS 2)
python_test(FILE long_range_actors.py MAX_NUM_PROC 1 GPU_SLOTS 1 SUFFIX 1_core)
python_test(FILE tabulated.py MAX_NUM_PROC 2)
python_test(FILE particle_slice.py MAX_NUM_PROC 4)
python_test(FILE rigid_bond.py MAX_NUM_PROC 4)
python_test(FILE rotation_per_particle.py MAX_NUM_PROC 4)
python_test(FILE rotational_inertia.py MAX_NUM_PROC 4)
python_test(FILE rotational-diffusion-aniso.py MAX_NUM_PROC 1 LABELS long)
python_test(FILE rotational_dynamics.py MAX_NUM_PROC 1)
python_test(FILE script_interface.py MAX_NUM_PROC 4)
python_test(FILE reaction_methods_interface.py MAX_NUM_PROC 2)
python_test(FILE reaction_ensemble.py MAX_NUM_PROC 4)
python_test(FILE reaction_complex.py MAX_NUM_PROC 1)
python_test(FILE reaction_bookkeeping.py MAX_NUM_PROC 1)
python_test(FILE widom_insertion.py MAX_NUM_PROC 1)
python_test(FILE constant_pH.py MAX_NUM_PROC 1)
python_test(FILE constant_pH_stats.py MAX_NUM_PROC 4 LABELS long)
python_test(FILE canonical_ensemble.py MAX_NUM_PROC 2)
python_test(FILE writevtf.py MAX_NUM_PROC 4)
# python_test(FILE lb_stokes_sphere.py MAX_NUM_PROC 4 GPU_SLOTS 1 LABELS long)
python_test(FILE lb_pressure_tensor.py MAX_NUM_PROC 1 GPU_SLOTS 3 LABELS long)
# python_test(FILE ek_fluctuations.py MAX_NUM_PROC 1) # TODO
python_test(FILE exclusions.py MAX_NUM_PROC 2)
python_test(FILE langevin_thermostat.py MAX_NUM_PROC 1)
python_test(FILE langevin_thermostat_stats.py MAX_NUM_PROC 1 LABELS long)
python_test(FILE brownian_dynamics.py MAX_NUM_PROC 1)
python_test(FILE brownian_dynamics_stats.py MAX_NUM_PROC 1 LABELS long)
python_test(FILE lees_edwards.py MAX_NUM_PROC 4)
python_test(FILE lb_lees_edwards.py MAX_NUM_PROC 1)
python_test(FILE lb_lees_edwards_particle_coupling.py MAX_NUM_PROC 1)
python_test(FILE lb_planar_couette.py MAX_NUM_PROC 1)
python_test(FILE nsquare.py MAX_NUM_PROC 4)
python_test(FILE virtual_sites_relative.py MAX_NUM_PROC 2)
python_test(FILE virtual_sites_relative_pbc.py MAX_NUM_PROC 2)
python_test(FILE virtual_sites_tracers.py MAX_NUM_PROC 2 DEPENDENCIES
            virtual_sites_tracers_common.py)
# python_test(FILE virtual_sites_tracers_gpu.py MAX_NUM_PROC 2 GPU_SLOTS 1
# DEPENDENCIES virtual_sites_tracers_common.py)
python_test(FILE regular_decomposition.py MAX_NUM_PROC 4)
python_test(FILE hybrid_decomposition.py MAX_NUM_PROC 1 SUFFIX 1_core)
python_test(FILE hybrid_decomposition.py MAX_NUM_PROC 4)
python_test(FILE integrator_npt.py MAX_NUM_PROC 4)
python_test(FILE integrator_npt_stats.py MAX_NUM_PROC 4 LABELS long)
python_test(FILE integrator_steepest_descent.py MAX_NUM_PROC 4)
python_test(FILE ibm.py MAX_NUM_PROC 2)
python_test(FILE dipolar_mdlc_p3m_scafacos_p2nfft.py MAX_NUM_PROC 1)
python_test(FILE dipolar_direct_summation.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE dipolar_p3m.py MAX_NUM_PROC 2)
python_test(FILE dipolar_interface.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE coulomb_interface.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb_stats.py MAX_NUM_PROC 2 GPU_SLOTS 2 LABELS long)
python_test(FILE lb_stats.py MAX_NUM_PROC 1 GPU_SLOTS 2 LABELS long SUFFIX
            1_core)
python_test(FILE force_cap.py MAX_NUM_PROC 2)
python_test(FILE dpd.py MAX_NUM_PROC 4)
python_test(FILE dpd_stats.py MAX_NUM_PROC 4 LABELS long)
python_test(FILE hat.py MAX_NUM_PROC 4)
python_test(FILE analyze_energy.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE analyze_mass_related.py MAX_NUM_PROC 4)
python_test(FILE rdf.py MAX_NUM_PROC 1)
python_test(FILE sf_simple_lattice.py MAX_NUM_PROC 1)
python_test(FILE coulomb_mixed_periodicity.py MAX_NUM_PROC 4)
python_test(FILE coulomb_cloud_wall_duplicated.py MAX_NUM_PROC 4 GPU_SLOTS 3)
python_test(FILE collision_detection.py MAX_NUM_PROC 4)
python_test(FILE collision_detection_interface.py MAX_NUM_PROC 2)
python_test(FILE lj.py MAX_NUM_PROC 4)
python_test(FILE dipole_field_tracking.py MAX_NUM_PROC 2)
python_test(FILE pairs.py MAX_NUM_PROC 4)
python_test(FILE polymer_linear.py MAX_NUM_PROC 4)
python_test(FILE polymer_diamond.py MAX_NUM_PROC 4)
python_test(FILE auto_exclusions.py MAX_NUM_PROC 4)
python_test(FILE observable_cylindrical.py MAX_NUM_PROC 4)
python_test(FILE observable_cylindricalLB.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE analyze_chains.py MAX_NUM_PROC 1)
python_test(FILE analyze_distance.py MAX_NUM_PROC 2)
python_test(FILE analyze_acf.py MAX_NUM_PROC 1)
python_test(FILE comfixed.py MAX_NUM_PROC 2)
python_test(FILE rescale.py MAX_NUM_PROC 2)
python_test(FILE array_properties.py MAX_NUM_PROC 4)
python_test(FILE analyze_distribution.py MAX_NUM_PROC 1)
python_test(FILE observable_profile.py MAX_NUM_PROC 4)
python_test(FILE observable_profileLB.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE rotate_system.py MAX_NUM_PROC 4)
python_test(FILE es_math.py MAX_NUM_PROC 1)
python_test(FILE random_pairs.py MAX_NUM_PROC 4)
python_test(FILE lb_electrohydrodynamics.py MAX_NUM_PROC 4 GPU_SLOTS 1)
python_test(FILE cluster_analysis.py MAX_NUM_PROC 4)
python_test(FILE pair_criteria.py MAX_NUM_PROC 4)
python_test(FILE actor.py MAX_NUM_PROC 1)
python_test(FILE drude.py MAX_NUM_PROC 2)
python_test(FILE thermostats_anisotropic.py MAX_NUM_PROC 4)
python_test(FILE thermalized_bond.py MAX_NUM_PROC 4)
python_test(FILE thole.py MAX_NUM_PROC 4)
python_test(FILE lb_slice.py MAX_NUM_PROC 2)
python_test(FILE lb_boundary_velocity.py MAX_NUM_PROC 1)
# python_test(FILE lb_boundary_volume_force.py MAX_NUM_PROC 2) # TODO
python_test(FILE lb_circular_couette.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb_thermo_virtual.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb_poiseuille.py MAX_NUM_PROC 4 GPU_SLOTS 1)
python_test(FILE lb_poiseuille_cylinder.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb_interpolation.py MAX_NUM_PROC 4 GPU_SLOTS 1)
python_test(FILE analyze_gyration_tensor.py MAX_NUM_PROC 1)
python_test(FILE oif_volume_conservation.py MAX_NUM_PROC 2)
python_test(FILE simple_pore.py MAX_NUM_PROC 1)
python_test(FILE field_test.py MAX_NUM_PROC 1)
python_test(FILE lb_boundary.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb_streaming.py MAX_NUM_PROC 4 GPU_SLOTS 1)
python_test(FILE lb_shear.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE lb_thermostat.py MAX_NUM_PROC 2 GPU_SLOTS 1)
# python_test(FILE lb_buoyancy_force.py MAX_NUM_PROC 2 GPU_SLOTS 1) # TODO
python_test(FILE lb_momentum_conservation.py MAX_NUM_PROC 2 GPU_SLOTS 1 LABELS
            long)
python_test(FILE lb_momentum_conservation.py MAX_NUM_PROC 1 GPU_SLOTS 1 LABELS
            long SUFFIX 1_core)
python_test(FILE lb_mass_conservation.py MAX_NUM_PROC 2)
python_test(FILE p3m_electrostatic_pressure.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE p3m_madelung.py MAX_NUM_PROC 2 GPU_SLOTS 2 LABELS long)
python_test(FILE sigint.py DEPENDENCIES sigint_child.py NO_MPI)
python_test(FILE observable_chain.py MAX_NUM_PROC 4)
python_test(FILE mpiio.py MAX_NUM_PROC 4)
python_test(FILE mpiio_exceptions.py MAX_NUM_PROC 1)
python_test(FILE gpu_availability.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE features.py MAX_NUM_PROC 1)
python_test(FILE decorators.py MAX_NUM_PROC 1)
python_test(FILE galilei.py MAX_NUM_PROC 4)
python_test(FILE linear_momentum.py MAX_NUM_PROC 4)
python_test(FILE linear_momentum_lb.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE mmm1d.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE stokesian_dynamics.py MAX_NUM_PROC 2 LABELS long)
python_test(FILE stokesian_thermostat.py MAX_NUM_PROC 2)
python_test(FILE elc.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE elc_vs_analytic.py MAX_NUM_PROC 2 GPU_SLOTS 1)
python_test(FILE rotation.py MAX_NUM_PROC 1)
python_test(FILE shapes.py MAX_NUM_PROC 1)
python_test(FILE h5md.py MAX_NUM_PROC 2)
python_test(FILE h5md.py MAX_NUM_PROC 1 SUFFIX 1_core)
python_test(FILE p3m_fft.py MAX_NUM_PROC 6)
if(${ESPRESSO_TEST_NP} GREATER_EQUAL 8)
  python_test(FILE p3m_fft.py MAX_NUM_PROC 8 SUFFIX 8_cores)
endif()
python_test(FILE p3m_tuning_exceptions.py MAX_NUM_PROC 1 GPU_SLOTS 1)
python_test(FILE integrator_exceptions.py MAX_NUM_PROC 1)
python_test(FILE utils.py MAX_NUM_PROC 1)
python_test(FILE npt_thermostat.py MAX_NUM_PROC 4)
python_test(FILE box_geometry.py MAX_NUM_PROC 1)
python_test(FILE lattice.py MAX_NUM_PROC 4)
python_test(FILE lattice_vtk.py MAX_NUM_PROC 4)
if(${ESPRESSO_TEST_NP} GREATER_EQUAL 6)
  python_test(FILE lattice_vtk.py MAX_NUM_PROC 6 SUFFIX 6_cores)
endif()
python_test(FILE ek_interface.py MAX_NUM_PROC 2)
python_test(FILE ek_diffusion.py MAX_NUM_PROC 1)
python_test(FILE ek_noflux.py MAX_NUM_PROC 1)
python_test(FILE ek_eof.py MAX_NUM_PROC 1)
python_test(FILE ek_fixedflux.py MAX_NUM_PROC 1)
python_test(FILE ek_bulk_reactions.py MAX_NUM_PROC 1)
python_test(FILE ek_indexed_reactions.py MAX_NUM_PROC 1)
python_test(FILE ek_fixeddensity.py MAX_NUM_PROC 1)
python_test(FILE ek_boundary.py MAX_NUM_PROC 2)
python_test(FILE ek_slice.py MAX_NUM_PROC 2)

set(ESPRESSO_CTEST_RESOURCE_SPEC_FILE resources.json)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${ESPRESSO_CTEST_RESOURCE_SPEC_FILE}
  ${CMAKE_CURRENT_BINARY_DIR}/${ESPRESSO_CTEST_RESOURCE_SPEC_FILE} COPYONLY)
add_custom_target(
  python_test_data
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
          ${CMAKE_CURRENT_BINARY_DIR}/data
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/tests_common.py
          ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/unittest_decorators.py
    ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/thermostats_common.py
    ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(
  check_python_parallel_odd
  COMMAND
    ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT} -L parallel_odd
    --resource-spec-file ${ESPRESSO_CTEST_RESOURCE_SPEC_FILE}
    ${ESPRESSO_CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_parallel_odd pypresso python_test_data)
add_custom_target(
  check_python_gpu
  COMMAND
    ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT} -L gpu
    --resource-spec-file ${ESPRESSO_CTEST_RESOURCE_SPEC_FILE}
    ${ESPRESSO_CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_gpu pypresso python_test_data)

add_custom_target(
  check_python_skip_long
  COMMAND
    ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT} -LE long
    --resource-spec-file ${ESPRESSO_CTEST_RESOURCE_SPEC_FILE}
    ${ESPRESSO_CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_skip_long pypresso python_test_data)

add_custom_target(
  check_python
  COMMAND
    ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
    --resource-spec-file ${ESPRESSO_CTEST_RESOURCE_SPEC_FILE}
    ${ESPRESSO_CTEST_ARGS} --output-on-failure)
add_dependencies(check_python pypresso python_test_data)
add_dependencies(check check_python)
