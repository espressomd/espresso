function(PYTHON_TEST)
  cmake_parse_arguments(TEST "" "FILE;MAX_NUM_PROC;SUFFIX"
                        "DEPENDS;DEPENDENCIES;LABELS" ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_${TEST_SUFFIX}")
  endif()
  set(TEST_FILE_CONFIGURED "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py")
  configure_file(${TEST_FILE} ${TEST_FILE_CONFIGURED})
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency})
  endforeach(dependency)
  set(TEST_FILE ${TEST_FILE_CONFIGURED})

  if(NOT DEFINED TEST_MAX_NUM_PROC)
    set(TEST_MAX_NUM_PROC ${TEST_NP})
  endif()

  if(${TEST_MAX_NUM_PROC} GREATER ${TEST_NP})
    set(TEST_NUM_PROC ${TEST_NP})
  else()
    set(TEST_NUM_PROC ${TEST_MAX_NUM_PROC})
  endif()

  if(EXISTS ${MPIEXEC})
    add_test(
      NAME ${TEST_NAME}
      COMMAND
        ${MPIEXEC} ${MPIEXEC_OVERSUBSCRIBE} ${MPIEXEC_NUMPROC_FLAG}
        ${TEST_NUM_PROC} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/pypresso
        ${PYPRESSO_OPTIONS} ${TEST_FILE} ${MPIEXEC_POSTFLAGS})
  else()
    add_test(${TEST_NAME} ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
             ${TEST_FILE})
  endif()
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${TEST_NUM_PROC}
                                               DEPENDS "${TEST_DEPENDS}")

  if("gpu" IN_LIST TEST_LABELS AND WITH_CUDA)
    set_tests_properties(${TEST_NAME} PROPERTIES RESOURCE_LOCK GPU)
  endif()

  if(${TEST_MAX_NUM_PROC} EQUAL 1)
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}")
  elseif(${TEST_MAX_NUM_PROC} EQUAL 2)
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS
                                                 "${TEST_LABELS};parallel")
  else()
    set_tests_properties(
      ${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS};parallel;parallel_odd")
  endif()

  set(python_tests ${python_tests} ${TEST_FILE} PARENT_SCOPE)
endfunction(PYTHON_TEST)

# Checkpointing tests: semicolon-separated list of mutually-compatible features.
# Separate features with hyphens, use a period to add an optional flag.
foreach(
  TEST_COMBINATION
  lb.active.walberla-p3m.cpu-lj-therm.lb;lb.active.walberla-p3m.elc-lj-therm.lb;lb.off-therm.npt-int.npt;lb.off-int.sd;lb.off-dp3m.cpu-therm.langevin-int.nvt;lb.off-therm.dpd-int.nvt;lb.off-scafacos-therm.bd-int.bd;lb.off-therm.sdm-int.sdm;lb.active.walberla-p3m.cpu-lj-therm.lb-mpi.1.core
)
  if(${TEST_COMBINATION} MATCHES "mpi\\.1\\.core")
    set(TEST_NPROCS 1)
  else()
    set(TEST_NPROCS 2)
  endif()
  if(${TEST_COMBINATION} MATCHES "\\.gpu")
    set(TEST_LABELS "gpu")
  else()
    set(TEST_LABELS "")
  endif()
  if(${TEST_COMBINATION} MATCHES "lb\\.off")
    set(TEST_BINARY_LIST "1")
  else()
    set(TEST_BINARY_LIST "1;0")
  endif()
  foreach(TEST_BINARY ${TEST_BINARY_LIST})
    python_test(FILE save_checkpoint.py MAX_NUM_PROC ${TEST_NPROCS} LABELS
                ${TEST_LABELS} SUFFIX ${TEST_COMBINATION}_${TEST_BINARY})
    python_test(
      FILE test_checkpoint.py MAX_NUM_PROC ${TEST_NPROCS} LABELS ${TEST_LABELS}
      SUFFIX ${TEST_COMBINATION}_${TEST_BINARY} DEPENDS
      save_checkpoint_${TEST_COMBINATION}_${TEST_BINARY})
  endforeach(TEST_BINARY)
endforeach(TEST_COMBINATION)
python_test(FILE cellsystem.py MAX_NUM_PROC 4)
python_test(FILE tune_skin.py MAX_NUM_PROC 1)
python_test(FILE constraint_homogeneous_magnetic_field.py MAX_NUM_PROC 4)
python_test(FILE constraint_shape_based.py MAX_NUM_PROC 2)
python_test(FILE coulomb_cloud_wall.py MAX_NUM_PROC 4 LABELS gpu)
python_test(FILE coulomb_tuning.py MAX_NUM_PROC 4 LABELS gpu long)
python_test(FILE accumulator_correlator.py MAX_NUM_PROC 4)
python_test(FILE accumulator_mean_variance.py MAX_NUM_PROC 4)
python_test(FILE accumulator_time_series.py MAX_NUM_PROC 1)
python_test(FILE dawaanr-and-dds-gpu.py MAX_NUM_PROC 1 LABELS gpu)
python_test(FILE dawaanr-and-bh-gpu.py MAX_NUM_PROC 1 LABELS gpu)
python_test(FILE dds-and-bh-gpu.py MAX_NUM_PROC 4 LABELS gpu)
python_test(FILE electrostaticInteractions.py MAX_NUM_PROC 2)
python_test(FILE engine_langevin.py MAX_NUM_PROC 4)
python_test(FILE engine_lb.py MAX_NUM_PROC 2)
python_test(FILE engine_lb_nsquare.py MAX_NUM_PROC 1)
python_test(FILE experimental_decorator.py)
python_test(FILE icc.py MAX_NUM_PROC 4)
python_test(FILE mass-and-rinertia_per_particle.py MAX_NUM_PROC 2 LABELS long)
python_test(FILE integrate.py MAX_NUM_PROC 4)
python_test(FILE interactions_bond_angle.py MAX_NUM_PROC 4)
python_test(FILE interactions_bonded_interface.py MAX_NUM_PROC 4)
python_test(FILE interactions_bonded.py MAX_NUM_PROC 2)
python_test(FILE interactions_dihedral.py MAX_NUM_PROC 4)
python_test(FILE interactions_non-bonded_interface.py MAX_NUM_PROC 4)
python_test(FILE interactions_non-bonded.py MAX_NUM_PROC 4)
python_test(FILE observables.py MAX_NUM_PROC 4)
python_test(FILE p3m_gpu.py MAX_NUM_PROC 4 LABELS gpu)
python_test(FILE particle.py MAX_NUM_PROC 4)
python_test(FILE pressure.py MAX_NUM_PROC 4)
python_test(FILE scafacos_dipoles_1d_2d.py MAX_NUM_PROC 4)
python_test(FILE scafacos_interface.py MAX_NUM_PROC 2)
python_test(FILE tabulated.py MAX_NUM_PROC 2)
python_test(FILE particle_slice.py MAX_NUM_PROC 4)
python_test(FILE rigid_bond.py MAX_NUM_PROC 4)
python_test(FILE rotation_per_particle.py MAX_NUM_PROC 4)
python_test(FILE rotational_inertia.py MAX_NUM_PROC 4)
python_test(FILE rotational-diffusion-aniso.py MAX_NUM_PROC 1 LABELS long)
python_test(FILE rotational_dynamics.py MAX_NUM_PROC 1)
python_test(FILE script_interface.py MAX_NUM_PROC 4)
python_test(FILE reaction_ensemble.py MAX_NUM_PROC 4)
python_test(FILE widom_insertion.py MAX_NUM_PROC 1)
python_test(FILE constant_pH.py MAX_NUM_PROC 1)
python_test(FILE constant_pH_stats.py MAX_NUM_PROC 4 LABELS long)
python_test(FILE writevtf.py MAX_NUM_PROC 4)
# python_test(FILE lb_stokes_sphere.py MAX_NUM_PROC 8 LABELS gpu long) # TODO WALBERLA
python_test(FILE lb_pressure_tensor.py MAX_NUM_PROC 1 LABELS gpu long)
# python_test(FILE ek_fluctuations.py MAX_NUM_PROC 1 LABELS gpu) # TODO WALBERLA
# python_test(FILE ek_charged_plate.py MAX_NUM_PROC 1 LABELS gpu) # TODO WALBERLA
# python_test(FILE ek_eof_one_species.py MAX_NUM_PROC 1 LABELS gpu SUFFIX x) # TODO WALBERLA
# python_test(FILE ek_eof_one_species.py MAX_NUM_PROC 1 LABELS gpu SUFFIX y) # TODO WALBERLA
# python_test(FILE ek_eof_one_species.py MAX_NUM_PROC 1 LABELS gpu SUFFIX z) # TODO WALBERLA
python_test(FILE exclusions.py MAX_NUM_PROC 2)
python_test(FILE langevin_thermostat.py MAX_NUM_PROC 1)
python_test(FILE langevin_thermostat_stats.py MAX_NUM_PROC 1 LABELS long)
python_test(FILE brownian_dynamics.py MAX_NUM_PROC 1)
python_test(FILE brownian_dynamics_stats.py MAX_NUM_PROC 1 LABELS long)
python_test(FILE lees_edwards.py MAX_NUM_PROC 1)
python_test(FILE lb_lees_edwards.py MAX_NUM_PROC 1)
python_test(FILE lb_lees_edwards_particle_coupling.py MAX_NUM_PROC 1)
python_test(FILE lb_planar_couette.py MAX_NUM_PROC 1)
python_test(FILE nsquare.py MAX_NUM_PROC 4)
python_test(FILE virtual_sites_relative.py MAX_NUM_PROC 2)
python_test(FILE virtual_sites_tracers_walberla.py MAX_NUM_PROC 1)
python_test(FILE domain_decomposition.py MAX_NUM_PROC 4)
# python_test(FILE integrator_npt.py MAX_NUM_PROC 4)  # TODO WALBERLA: Random segfault in the CellSystem
# python_test(FILE integrator_npt_stats.py MAX_NUM_PROC 4 LABELS long)  # TODO WALBERLA: Random segfault in the CellSystem
python_test(FILE integrator_steepest_descent.py MAX_NUM_PROC 4)
python_test(FILE ibm.py MAX_NUM_PROC 2)
python_test(FILE dipolar_mdlc_p3m_scafacos_p2nfft.py MAX_NUM_PROC 1)
python_test(FILE dipolar_direct_summation.py MAX_NUM_PROC 1 LABELS gpu)
python_test(FILE dipolar_p3m.py MAX_NUM_PROC 2)
python_test(FILE dipolar_interface.py MAX_NUM_PROC 1 LABELS gpu)
python_test(FILE dipolar_mpi_exceptions.py MAX_NUM_PROC 2)
python_test(FILE lattice.py MAX_NUM_PROC 4)
python_test(FILE lb.py MAX_NUM_PROC 2)
python_test(FILE lb_walberla.py MAX_NUM_PROC 2)
python_test(FILE lb_stats.py MAX_NUM_PROC 2)
python_test(FILE lb_vtk.py MAX_NUM_PROC 1)
python_test(FILE force_cap.py MAX_NUM_PROC 2)
python_test(FILE dpd.py MAX_NUM_PROC 4)
python_test(FILE dpd_stats.py MAX_NUM_PROC 4 LABELS long)
python_test(FILE hat.py MAX_NUM_PROC 4)
python_test(FILE analyze_energy.py MAX_NUM_PROC 2)
python_test(FILE analyze_mass_related.py MAX_NUM_PROC 4)
python_test(FILE rdf.py MAX_NUM_PROC 1)
python_test(FILE sf_simple_lattice.py MAX_NUM_PROC 1)
python_test(FILE coulomb_mixed_periodicity.py MAX_NUM_PROC 4)
python_test(FILE coulomb_cloud_wall_duplicated.py MAX_NUM_PROC 4 LABELS gpu)
python_test(FILE collision_detection.py MAX_NUM_PROC 4)
python_test(FILE lj.py MAX_NUM_PROC 4)
python_test(FILE pairs.py MAX_NUM_PROC 4)
python_test(FILE polymer_linear.py MAX_NUM_PROC 4)
python_test(FILE polymer_diamond.py MAX_NUM_PROC 4)
python_test(FILE auto_exclusions.py MAX_NUM_PROC 1)
python_test(FILE observable_cylindrical.py MAX_NUM_PROC 4)
python_test(FILE observable_cylindricalLB.py MAX_NUM_PROC 1)
python_test(FILE analyze_chains.py MAX_NUM_PROC 1)
python_test(FILE analyze_distance.py MAX_NUM_PROC 1)
python_test(FILE analyze_acf.py MAX_NUM_PROC 1)
python_test(FILE comfixed.py MAX_NUM_PROC 2)
python_test(FILE rescale.py MAX_NUM_PROC 2)
python_test(FILE array_properties.py MAX_NUM_PROC 4)
python_test(FILE analyze_distribution.py MAX_NUM_PROC 1)
python_test(FILE observable_profile.py MAX_NUM_PROC 4)
python_test(FILE observable_profileLB.py MAX_NUM_PROC 1)
python_test(FILE rotate_system.py MAX_NUM_PROC 4)
python_test(FILE es_math.py MAX_NUM_PROC 1)
python_test(FILE random_pairs.py MAX_NUM_PROC 4)
python_test(FILE lb_electrohydrodynamics.py MAX_NUM_PROC 2)
python_test(FILE cluster_analysis.py MAX_NUM_PROC 4)
python_test(FILE pair_criteria.py MAX_NUM_PROC 4)
python_test(FILE actor.py MAX_NUM_PROC 1)
python_test(FILE drude.py MAX_NUM_PROC 2)
python_test(FILE thermalized_bond.py MAX_NUM_PROC 4)
python_test(FILE thole.py MAX_NUM_PROC 4)
python_test(FILE lb_slice.py MAX_NUM_PROC 1)
python_test(FILE lb_boundary_velocity.py MAX_NUM_PROC 1)
# python_test(FILE lb_boundary_volume_force.py MAX_NUM_PROC 2) # TODO WALBERLA
python_test(FILE lb_circular_couette.py MAX_NUM_PROC 2)
python_test(FILE lb_thermo_virtual.py MAX_NUM_PROC 2)
python_test(FILE lb_poiseuille.py MAX_NUM_PROC 2)
python_test(FILE lb_poiseuille_cylinder.py MAX_NUM_PROC 2)
# python_test(FILE lb_interpolation.py MAX_NUM_PROC 2) # TODO WALBERLA (#4381)
python_test(FILE analyze_gyration_tensor.py MAX_NUM_PROC 1)
python_test(FILE oif_volume_conservation.py MAX_NUM_PROC 2)
python_test(FILE simple_pore.py MAX_NUM_PROC 1)
python_test(FILE field_test.py MAX_NUM_PROC 1)
python_test(FILE lb_boundary.py MAX_NUM_PROC 2)
python_test(FILE lb_streaming.py MAX_NUM_PROC 2)
# python_test(FILE lb_shear.py MAX_NUM_PROC 2) # TODO WALBERLA (#4381)
python_test(FILE lb_thermostat.py MAX_NUM_PROC 1)
# python_test(FILE lb_buoyancy_force.py MAX_NUM_PROC 2) # TODO WALBERLA (dynamic UBB)
python_test(FILE lb_momentum_conservation.py MAX_NUM_PROC 2)
python_test(FILE lb_momentum_conservation_nsquare.py MAX_NUM_PROC 1)
python_test(FILE p3m_electrostatic_pressure.py MAX_NUM_PROC 2)
python_test(FILE sigint.py DEPENDENCIES sigint_child.py MAX_NUM_PROC 1)
python_test(FILE lb_density.py MAX_NUM_PROC 2)
python_test(FILE observable_chain.py MAX_NUM_PROC 4)
python_test(FILE mpiio.py MAX_NUM_PROC 4)
python_test(FILE gpu_availability.py MAX_NUM_PROC 2 LABELS gpu)
python_test(FILE features.py MAX_NUM_PROC 1)
python_test(FILE decorators.py MAX_NUM_PROC 1)
python_test(FILE galilei.py MAX_NUM_PROC 32)
python_test(FILE linear_momentum.py MAX_NUM_PROC 4)
python_test(FILE linear_momentum_lb.py MAX_NUM_PROC 2)
python_test(FILE mmm1d.py MAX_NUM_PROC 2)
python_test(FILE mmm1d_gpu.py MAX_NUM_PROC 1 LABELS gpu)
python_test(FILE stokesian_dynamics.py MAX_NUM_PROC 2 LABELS long)
python_test(FILE stokesian_thermostat.py MAX_NUM_PROC 2)
python_test(FILE elc.py MAX_NUM_PROC 2)
python_test(FILE elc_vs_analytic.py MAX_NUM_PROC 2)
python_test(FILE rotation.py MAX_NUM_PROC 1)
python_test(FILE shapes.py MAX_NUM_PROC 1)
python_test(FILE h5md.py MAX_NUM_PROC 2)
python_test(FILE mdanalysis.py MAX_NUM_PROC 2)
python_test(FILE p3m_fft.py MAX_NUM_PROC 6)
if(${TEST_NP} GREATER_EQUAL 8)
  python_test(FILE p3m_fft.py MAX_NUM_PROC 8 SUFFIX 8_cores)
endif()
python_test(FILE p3m_tuning_exceptions.py MAX_NUM_PROC 1 LABELS gpu)
python_test(FILE integrator_exceptions.py MAX_NUM_PROC 1)
python_test(FILE utils.py MAX_NUM_PROC 1)
# python_test(FILE npt_thermostat.py MAX_NUM_PROC 4) # TODO WALBERLA: Random segfault in the CellSystem
python_test(FILE box_geometry.py MAX_NUM_PROC 1)

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
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/virtual_sites_tracers_common.py
    ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/ek_common.py
          ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/thermostats_common.py
    ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(
  check_python_parallel_odd
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -L parallel_odd
          ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_parallel_odd pypresso python_test_data)
add_custom_target(
  check_python_gpu COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -L
                           gpu ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_gpu pypresso python_test_data)

add_custom_target(
  check_python_skip_long
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -LE long
          ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_skip_long pypresso python_test_data)

add_custom_target(
  check_python COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                       ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python pypresso python_test_data)
add_dependencies(check check_python)
