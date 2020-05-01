set(EspressoCore_SRC
    accumulators.cpp
    bond_error.cpp
    cells.cpp
    collision.cpp
    comfixed_global.cpp
    communication.cpp
    constraints.cpp
    debug.cpp
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
    npt.cpp
    nsquare.cpp
    partCfg_global.cpp
    particle_data.cpp
    polymer.cpp
    pressure.cpp
    rattle.cpp
    reaction_ensemble.cpp
    rotate_system.cpp
    rotation.cpp
    RuntimeErrorCollector.cpp
    RuntimeError.cpp
    RuntimeErrorStream.cpp
    statistics_chain.cpp
    statistics.cpp
    SystemInterface.cpp
    thermostat.cpp
    tuning.cpp
    virtual_sites.cpp
    exclusions.cpp
    CellStructure.cpp
    PartCfg.cpp
    AtomDecomposition.cpp)

if(CUDA)
  set(EspressoCuda_SRC
      actor/DipolarBarnesHut_cuda.cu
      actor/DipolarDirectSum_cuda.cu
      actor/Mmm1dgpuForce_cuda.cu
      cuda_common_cuda.cu
      cuda_init.cpp
      cuda_init_cuda.cu
      cuda_interface.cpp
      CudaHostAllocator.cu
      electrostatics_magnetostatics/p3m_gpu_cuda.cu
      electrostatics_magnetostatics/p3m_gpu_error_cuda.cu
      EspressoSystemInterface_cuda.cu
      grid_based_algorithms/electrokinetics_cuda.cu
      grid_based_algorithms/lbgpu_cuda.cu
      grid_based_algorithms/fd-electrostatics_cuda.cu
      grid_based_algorithms/electrokinetics.cpp
      grid_based_algorithms/lbgpu.cpp
      virtual_sites/lb_inertialess_tracers_cuda.cu)

  add_gpu_library(EspressoCore SHARED ${EspressoCore_SRC} ${EspressoCuda_SRC})
else(CUDA)
  add_library(EspressoCore SHARED ${EspressoCore_SRC})
endif(CUDA)

install(TARGETS EspressoCore LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

target_link_libraries(
  EspressoCore
  PRIVATE EspressoConfig shapes Profiler
          "$<$<BOOL:${FFTW3_FOUND}>:${FFTW3_LIBRARIES}>"
          $<$<BOOL:${SCAFACOS}>:Scafacos> cxx_interface
  PUBLIC utils MPI::MPI_CXX Random123 Boost::serialization Boost::mpi
         "$<$<BOOL:${H5MD}>:${HDF5_LIBRARIES}>"
         $<$<BOOL:${H5MD}>:Boost::filesystem> $<$<BOOL:${H5MD}>:h5xx>)

target_include_directories(
  EspressoCore
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
         "$<$<BOOL:${H5MD}>:${CMAKE_CURRRENT_SOURCE_DIR}/io/writer>" SYSTEM
  PUBLIC "$<$<BOOL:${H5MD}>:${HDF5_INCLUDE_DIRS}>" SYSTEM
  PRIVATE "$<$<BOOL:${FFTW3_FOUND}>:${FFTW3_INCLUDE_DIR}>")

target_compile_definitions(EspressoCore PUBLIC $<$<BOOL:${H5MD}>:H5XX_USE_MPI>)

add_subdirectory(accumulators)
add_subdirectory(actor)
add_subdirectory(bonded_interactions)
add_subdirectory(cluster_analysis)
add_subdirectory(constraints)
add_subdirectory(electrostatics_magnetostatics)
add_subdirectory(grid_based_algorithms)
add_subdirectory(immersed_boundary)
add_subdirectory(integrators)
add_subdirectory(io)
add_subdirectory(nonbonded_interactions)
add_subdirectory(object-in-fluid)
add_subdirectory(observables)
add_subdirectory(virtual_sites)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)
