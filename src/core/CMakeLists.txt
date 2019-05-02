set(EspressoCore_SRC
  accumulators.cpp
  cells.cpp
  collision.cpp
  comfixed_global.cpp
  communication.cpp
  constraints.cpp
  debug.cpp
  diamond.cpp
  domain_decomposition.cpp
  dpd.cpp
  energy.cpp
  errorhandling.cpp
  EspressoSystemInterface.cpp
  forcecap.cpp
  forces.cpp
  galilei.cpp
  ghosts.cpp
  global.cpp
  grid.cpp
  immersed_boundaries.cpp
  event.cpp
  integrate.cpp
  layered.cpp
  metadynamics.cpp
  minimize_energy.cpp
  npt.cpp
  nsquare.cpp
  partCfg_global.cpp
  particle_data.cpp
  polymer.cpp
  polynom.cpp
  pressure.cpp
  random.cpp
  rattle.cpp
  reaction_ensemble.cpp
  rotate_system.cpp
  rotation.cpp
  RuntimeErrorCollector.cpp
  RuntimeError.cpp
  RuntimeErrorStream.cpp
  specfunc.cpp
  statistics_chain.cpp
  statistics.cpp
  swimmer_reaction.cpp
  SystemInterface.cpp
  thermostat.cpp
  tuning.cpp
  virtual_sites.cpp
  accumulators/Correlator.cpp
  accumulators/MeanVarianceCalculator.cpp
  actor/ActorList.cpp
  actor/DipolarBarnesHut.cpp
  actor/DipolarDirectSum.cpp
  actor/HarmonicOrientationWell.cpp
  actor/HarmonicWell.cpp
  actor/Mmm1dgpuForce.cpp
  bonded_interactions/angle_cosine.cpp
  bonded_interactions/angle_cossquare.cpp
  bonded_interactions/angle_harmonic.cpp
  bonded_interactions/bonded_coulomb.cpp
  bonded_interactions/bonded_coulomb_sr.cpp
  bonded_interactions/bonded_interaction_data.cpp
  bonded_interactions/bonded_tab.cpp
  bonded_interactions/dihedral.cpp
  bonded_interactions/fene.cpp
  bonded_interactions/harmonic.cpp
  bonded_interactions/harmonic_dumbbell.cpp
  bonded_interactions/quartic.cpp
  bonded_interactions/subt_lj.cpp
  bonded_interactions/thermalized_bond.cpp
  bonded_interactions/umbrella.cpp
  constraints/HomogeneousMagneticField.cpp
  constraints/ShapeBasedConstraint.cpp
  electrostatics_magnetostatics/debye_hueckel.cpp
  electrostatics_magnetostatics/elc.cpp
  electrostatics_magnetostatics/icc.cpp
  electrostatics_magnetostatics/magnetic_non_p3m_methods.cpp
  electrostatics_magnetostatics/mdlc_correction.cpp
  electrostatics_magnetostatics/mmm1d.cpp
  electrostatics_magnetostatics/mmm2d.cpp
  electrostatics_magnetostatics/mmm-common.cpp
  electrostatics_magnetostatics/p3m-common.cpp
  electrostatics_magnetostatics/p3m.cpp
  electrostatics_magnetostatics/p3m-dipolar.cpp
  electrostatics_magnetostatics/p3m_gpu.cpp
  electrostatics_magnetostatics/scafacos.cpp
  electrostatics_magnetostatics/fft.cpp
  electrostatics_magnetostatics/coulomb.cpp
  electrostatics_magnetostatics/dipole.cpp
  electrostatics_magnetostatics/reaction_field.cpp
  grid_based_algorithms/halo.cpp
  grid_based_algorithms/lattice.cpp
  grid_based_algorithms/lb_boundaries.cpp
  grid_based_algorithms/lb.cpp
  grid_based_algorithms/lb_interface.cpp
  grid_based_algorithms/lb_interpolation.cpp
  grid_based_algorithms/lb_particle_coupling.cpp
  immersed_boundary/ibm_tribend.cpp
  immersed_boundary/ibm_triel.cpp
  immersed_boundary/ImmersedBoundaries.cpp
  nonbonded_interactions/bmhtf-nacl.cpp
  nonbonded_interactions/buckingham.cpp
  nonbonded_interactions/gaussian.cpp
  nonbonded_interactions/gb.cpp
  nonbonded_interactions/hat.cpp
  nonbonded_interactions/hertzian.cpp
  nonbonded_interactions/ljcos2.cpp
  nonbonded_interactions/ljcos.cpp
  nonbonded_interactions/lj.cpp
  nonbonded_interactions/ljgen.cpp
  nonbonded_interactions/morse.cpp
  nonbonded_interactions/nonbonded_interaction_data.cpp
  nonbonded_interactions/nonbonded_tab.cpp
  nonbonded_interactions/soft_sphere.cpp
  nonbonded_interactions/steppot.cpp
  nonbonded_interactions/thole.cpp
  nonbonded_interactions/wca.cpp
  object-in-fluid/affinity.cpp
  object-in-fluid/membrane_collision.cpp
  object-in-fluid/oif_global_forces.cpp
  object-in-fluid/oif_local_forces.cpp
  object-in-fluid/out_direction.cpp
  observables/CylindricalLBFluxDensityProfileAtParticlePositions.cpp
  observables/CylindricalLBVelocityProfileAtParticlePositions.cpp
  observables/CylindricalLBVelocityProfile.cpp
  observables/LBVelocityProfile.cpp
  virtual_sites/lb_inertialess_tracers.cpp
  virtual_sites/lb_inertialess_tracers_cuda_interface.cpp
  virtual_sites/virtual_sites_com.cpp
  virtual_sites/VirtualSitesInertialessTracers.cpp
  virtual_sites/VirtualSitesRelative.cpp
  )

if(CUDA)
  set(EspressoCuda_SRC
    cuda_common_cuda.cu
    cuda_init.cpp
    cuda_init_cuda.cu
    cuda_interface.cpp
    EspressoSystemInterface_cuda.cu
    actor/DipolarBarnesHut_cuda.cu
    actor/DipolarDirectSum_cuda.cu
    actor/HarmonicOrientationWell_cuda.cu
    actor/HarmonicWell_cuda.cu
    actor/Mmm1dgpuForce_cuda.cu
    electrostatics_magnetostatics/p3m_gpu_cuda.cu
    electrostatics_magnetostatics/p3m_gpu_error_cuda.cu
    grid_based_algorithms/electrokinetics_cuda.cu
    grid_based_algorithms/lbgpu_cuda.cu
    grid_based_algorithms/fd-electrostatics_cuda.cu
    virtual_sites/lb_inertialess_tracers_cuda.cu
    grid_based_algorithms/electrokinetics.cpp
    grid_based_algorithms/lbgpu.cpp
    )

  add_gpu_library(EspressoCore SHARED ${EspressoCore_SRC} ${EspressoCuda_SRC})
else(CUDA)
  add_library(EspressoCore SHARED ${EspressoCore_SRC})
endif(CUDA)

install(TARGETS EspressoCore LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(EspressoCore PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(EspressoCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(WITH_COVERAGE)
  target_compile_options(EspressoCore PUBLIC "$<$<CONFIG:Release>:-g>")
  target_compile_options(EspressoCore PUBLIC "$<$<CONFIG:Release>:-O0>")
  target_compile_options(EspressoCore PUBLIC "$<$<CXX_COMPILER_ID:Clang>:-fprofile-instr-generate>")
  target_compile_options(EspressoCore PUBLIC "$<$<CXX_COMPILER_ID:Clang>:-fcoverage-mapping>")
  target_compile_options(EspressoCore PUBLIC "$<$<NOT:$<CXX_COMPILER_ID:Clang>>:--coverage>")
  target_compile_options(EspressoCore PUBLIC "$<$<NOT:$<CXX_COMPILER_ID:Clang>>:-fprofile-arcs>")
  target_compile_options(EspressoCore PUBLIC "$<$<NOT:$<CXX_COMPILER_ID:Clang>>:-ftest-coverage>")
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_link_libraries(EspressoCore PUBLIC gcov)
  endif()
endif()

target_link_libraries(EspressoCore PUBLIC EspressoConfig utils)
target_link_libraries(EspressoCore PUBLIC Boost::serialization Boost::mpi MPI::MPI_CXX Random123)
target_link_libraries(EspressoCore PRIVATE Profiler)
target_link_libraries(EspressoCore PRIVATE pdbparser)

if(FFTW3_FOUND)
  target_include_directories(EspressoCore PRIVATE SYSTEM ${FFTW3_INCLUDE_DIR})
  target_link_libraries(EspressoCore PRIVATE ${FFTW3_LIBRARIES})
endif()

if(SCAFACOS)
  target_link_libraries(EspressoCore PRIVATE Scafacos)
endif(SCAFACOS)

# Subdirectories
add_subdirectory(io)


if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)

add_subdirectory(grid_based_algorithms)
add_subdirectory(shapes)
add_subdirectory(cluster_analysis)
