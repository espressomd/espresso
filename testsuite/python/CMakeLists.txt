set(py_tests  bondedInteractions.py
              cellsystem.py
              constraint_homogeneous_magnetic_field.py
              constraint_shape_based.py
              coulomb_cloud_wall.py
              coulomb_cloud_wall_tune.py
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
              virtual_sites_relative.py
              script_interface_object_params.py
              lbgpu_remove_total_momentum.py
              tabulated.py
              reaction_ensemble.py
              constant_pH.py
              writevtf.py
              lb_stokes_sphere_gpu.py
              ek_eof_one_species_x.py
              ek_eof_one_species_y.py
              ek_eof_one_species_z.py
              ek_eof_one_species_x_nonlinear.py
              ek_eof_one_species_y_nonlinear.py
              ek_eof_one_species_z_nonlinear.py
              exclusions.py
)
if(PY_H5PY)
  set(py_tests ${py_tests} h5md.py)
endif(PY_H5PY)
add_custom_target(python_tests
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

foreach(testfile ${py_tests})
  get_filename_component(basename ${testfile} NAME_WE)
  add_test(${basename} ${CMAKE_BINARY_DIR}/pypresso ${testfile})
endforeach(testfile ${py_tests})

add_custom_target(check_python COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check_python pypresso python_tests)
add_dependencies(check check_python)
