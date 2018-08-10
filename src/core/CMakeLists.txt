set(EspressoCore_SRC
    accumulators.cpp
    cells.cpp
    collision.cpp
    comfixed_global.cpp
    communication.cpp
    constraints.cpp
    cuda_init.cpp
    cuda_interface.cpp
    debug.cpp
    domain_decomposition.cpp
    dpd.cpp
    elc.cpp
    electrokinetics.cpp
    electrokinetics_pdb_parse.cpp
    energy.cpp
    errorhandling.cpp
    EspressoSystemInterface.cpp
    fft-common.cpp
    fft.cpp
    fft-dipolar.cpp
    forcecap.cpp
    forces.cpp
    galilei.cpp
    Geometry.cpp
    ghmc.cpp
    ghosts.cpp
    global.cpp
    grid.cpp
    halo.cpp
    iccp3m.cpp
    immersed_boundaries.cpp
    initialize.cpp
    integrate.cpp
    lattice.cpp
    layered.cpp
    lbboundaries.cpp
    lb.cpp
    lb-d3q19.cpp
    lbgpu.cpp
    maggs.cpp
    magnetic_non_p3m_methods.cpp
    mdlc_correction.cpp
    metadynamics.cpp
    minimize_energy.cpp
    mmm1d.cpp
    mmm2d.cpp
    mmm-common.cpp
    molforces.cpp
    MpiCallbacks.cpp
    nemd.cpp
    npt.cpp
    nsquare.cpp
    p3m-common.cpp
    p3m.cpp
    p3m-dipolar.cpp
    p3m_gpu.cpp
    partCfg_global.cpp
    particle_data.cpp
    PdbParser.cpp
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
    scafacos.cpp
    specfunc.cpp
    statistics_chain.cpp
    statistics_cluster.cpp
    statistics.cpp
    statistics_fluid.cpp
    swimmer_reaction.cpp
    SystemInterface.cpp
    thermostat.cpp
    topology.cpp
    tuning.cpp
    utils.cpp
    Vector.cpp
    virtual_sites.cpp
)


add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-features.hpp ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  ${CMAKE_SOURCE_DIR}/src/features.def
  ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
  ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  )
add_custom_target(
  generate_config_features
  DEPENDS ${CMAKE_BINARY_DIR}/src/core/config-features.hpp ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
)

add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
  COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)


add_library(EspressoConfig SHARED config-features.cpp)
add_dependencies(EspressoConfig myconfig check_myconfig generate_config_features)
install(TARGETS EspressoConfig LIBRARY DESTINATION ${LIBDIR})
set_target_properties(EspressoConfig PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(EspressoConfig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(EspressoCore SHARED ${EspressoCore_SRC} config-version.cpp)
add_dependencies(EspressoCore EspressoConfig)
install(TARGETS EspressoCore LIBRARY DESTINATION ${LIBDIR})
set_target_properties(EspressoCore PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(EspressoCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(EspressoCore PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
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

target_link_libraries(EspressoCore PRIVATE EspressoConfig ${LIBRARIES} Actor ObjectInFluid ImmersedBoundary Shapes Constraints Observables Accumulators cluster_analysis VirtualSites)


if(SCAFACOS)
  target_link_libraries(EspressoCore PRIVATE Scafacos)
endif(SCAFACOS)

if(GSL)
  target_include_directories(EspressoCore PUBLIC ${GSL_INCLUDE_DIRS})
endif(GSL)


if(GSL)
  target_link_libraries(EspressoCore PRIVATE gsl gslcblas)
endif(GSL)

if(CUDA)
    set(EspressoCuda_SRC
        cuda_common_cuda.cu
        cuda_init_cuda.cu
        electrokinetics_cuda.cu
        EspressoSystemInterface_cuda.cu
        fd-electrostatics_cuda.cu
        lbgpu_cuda.cu
        p3m_gpu_cuda.cu
        p3m_gpu_error_cuda.cu
    )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC})
  install(TARGETS EspressoCuda DESTINATION ${LIBDIR})
  set_target_properties(EspressoCuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(EspressoCuda EspressoConfig)

  CUDA_ADD_CUFFT_TO_TARGET(EspressoCuda)
  target_link_libraries(EspressoCore PRIVATE EspressoCuda)
endif(CUDA)

# Subdirectories
add_subdirectory(observables)
add_subdirectory(accumulators)
add_subdirectory(io)
if(HDF5_FOUND)
    target_link_libraries(EspressoCore PUBLIC H5mdCore)
endif(HDF5_FOUND)
target_link_libraries(EspressoCore PRIVATE mpiio)
target_link_libraries(EspressoCore PRIVATE pdbreader)
target_link_libraries(EspressoCore PRIVATE bonded_interactions)
target_link_libraries(EspressoCore PRIVATE nonbonded_interactions)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)

add_subdirectory(actor)
add_subdirectory(bonded_interactions)
add_subdirectory(immersed_boundary)
add_subdirectory(virtual_sites)
add_subdirectory(lbboundaries)
add_subdirectory(nonbonded_interactions)
add_subdirectory(object-in-fluid)
add_subdirectory(shapes)
add_subdirectory(constraints)
add_subdirectory(cluster_analysis)

if(SCAFACOS)
  add_subdirectory(scafacos)
endif(SCAFACOS)
