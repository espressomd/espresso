set(py_tests  bondedInteractions.py
              cellsystem.py
              constraint_homogeneous_magnetic_field.py
              constraint_shape_based.py
              coulomb_cloud_wall.py
              correlation.py
              dawaanr-and-dds-gpu.py
              electrostaticInteractions.py
              engine_langevin.py
              engine_lb.py
              engine_lbgpu.py
              ewald_gpu.py
              icc.py
              magnetostaticInteractions.py
              mass-and-rinertia_per_particle.py
              nonBondedInteractions.py
              observables.py
              p3m_gpu.py
              particle.py
              scafacos_dipoles_1d_2d.py
	            tabulated.py
              particle_slice.py
              rotational_inertia.py
              script_interface_object_params.py
              lbgpu_remove_total_momentum.py
              tabulated.py
              writevtf.py
              lb_stokes_sphere_gpu.py
              ek_eof_one_species_x.py
              ek_eof_one_species_y.py
              ek_eof_one_species_z.py
              ek_eof_one_species_x_nonlinear.py
              ek_eof_one_species_y_nonlinear.py
              ek_eof_one_species_z_nonlinear.py
)
if(PY_H5PY)
  set(py_tests ${py_tests} h5md.py)
endif(PY_H5PY)
add_custom_target(python_tests
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

if(NOT DEFINED TEST_NP)
  include(ProcessorCount)
  ProcessorCount(NP)
  math(EXPR TEST_NP "${NP}/2 + 1")
endif()

message("Np ${NP}")

foreach(testfile ${py_tests})
  get_filename_component(basename ${testfile} NAME_WE)

  if(EXISTS ${MPIEXEC})
    add_test(${basename} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${CMAKE_BINARY_DIR}/pypresso ${testfile})
  else()
    add_test(${basename} ${CMAKE_BINARY_DIR}/pypresso ${testfile})
  endif()
endforeach(testfile ${py_tests})

add_custom_target(check_python COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check_python pypresso python_tests)
add_dependencies(check check_python)
