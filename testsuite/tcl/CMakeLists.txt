set(tcl_tests  analysis.tcl
               angle.tcl
               bonded_coulomb.tcl
               collision-detection-angular.tcl
               collision-detection-centers.tcl
               collision-detection-glue.tcl
               collision-detection-poc.tcl
               comforce.tcl
               comfixed.tcl
               command_syntax.tcl
               constraints.tcl
               constraints_reflecting.tcl 
               correlation.tcl 
               constraints_rhomboid.tcl 
               coulomb_cloud_wall.tcl 
               dawaanr-and-dds-gpu.tcl 
               dh.tcl dielectric.tcl 
               dihedral.tcl dpd.tcl 
               ek_eof_one_species_x.tcl 
               ek_eof_one_species_y.tcl 
               ek_eof_one_species_z.tcl 
               ek_eof_one_species_x_nonlinear.tcl 
               ek_eof_one_species_y_nonlinear.tcl 
               ek_eof_one_species_z_nonlinear.tcl 
               ek_electrostatics_coupling_x.tcl 
               ek_electrostatics_coupling_y.tcl 
               ek_electrostatics_coupling_z.tcl 
               ek_reaction.tcl 
               el2d.tcl 
               el2d_nonneutral.tcl 
               engine_langevin.tcl 
               engine_lb.tcl 
               engine_lbgpu.tcl 
               exclusions.tcl 
               external_potential.tcl 
               fene.tcl 
               gb.tcl 
               ghmc.tcl 
               harm.tcl 
               quartic.tcl 
               iccp3m.tcl 
               immersed_boundary.tcl 
               immersed_boundary_gpu.tcl 
               intpbc.tcl 
               intppbc.tcl 
               kinetic.tcl 
               langevin.tcl 
               langevin_per_particle.tcl 
               layered.tcl 
               lb.tcl 
               lb_fluid_coupling.tcl 
               lb_fluid_coupling_gpu.tcl 
               lb_gpu.tcl 
               lb_planar.tcl 
               lb_planar_gpu.tcl 
               lb_planar_embedded_particles.tcl 
               lb_planar_embedded_particles_gpu.tcl 
               lb_stokes_sphere.tcl 
               lb_stokes_sphere_gpu.tcl 
               lees_edwards.tcl lj.tcl 
               lj-cos.tcl 
               lj-generic.tcl 
               madelung.tcl 
               maggs.tcl 
               magnetic-field.tcl 
               mass.tcl 
               mass-and-rinertia.tcl 
               mass-and-rinertia_per_particle.tcl 
               mdlc.tcl 
               minimize_energy.tcl 
               minimize_energy_rotation.tcl 
               mmm1d.tcl 
               mmm1dgpu.tcl 
               ewaldgpu.tcl 
               npt.tcl 
               nsquare.tcl 
               nve_pe.tcl 
               object_in_fluid.tcl 
               object_in_fluid_gpu.tcl 
               observable.tcl p3m.tcl 
               p3m_gpu.tcl 
               p3m_gpu_simple_noncubic.tcl 
               p3m_magnetostatics.tcl 
               p3m_magnetostatics2.tcl 
               p3m_simple_noncubic.tcl 
               pdb_parser.tcl 
               rotate-system.tcl 
               rotate-system-dipoles.tcl 
               rotation.tcl 
               rotation_per_particle.tcl 
               sc_fluid_gpu.tcl 
               sd_ewald.tcl 
               sd_two_spheres.tcl 
               sd_thermalization.tcl 
               tabulated.tcl 
               tunable_slip.tcl 
               uwerr.tcl 
               virtual-sites.tcl 
               virtual-sites-rotation.tcl)

add_custom_target(tcl_tests
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


if(WITH_MPI)
    include(ProcessorCount)
    ProcessorCount(NP)
    math(EXPR NP "${NP}/2 + 1")
endif(WITH_MPI)

if(WITH_TCL)
  foreach(testfile ${tcl_tests})
    get_filename_component(basename ${testfile} NAME_WE)
    set(TEST_EXECUTABLE ${PYTHON_EXECUTABLE} test_wrapper.py ${CMAKE_BINARY_DIR}/Espresso)
    if(WITH_MPI AND (EXISTS ${MPIEXEC}))
      add_test(${basename} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${TEST_EXECUTABLE} ${testfile})
    else()
      add_test(${basename} ${TEST_EXECUTABLE} ${testfile})
    endif()
    set_tests_properties(${basename} PROPERTIES SKIP_RETURN_CODE 42)
  endforeach(testfile ${tcl_tests})
  add_custom_target(check_tcl COMMAND ${CMAKE_CTEST_COMMAND}
    DEPENDS Espresso tcl_tests)
  add_dependencies(check check_tcl)
endif(WITH_TCL)

