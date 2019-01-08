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
    energy.cpp
    errorhandling.cpp
    EspressoSystemInterface.cpp
    fft-common.cpp
    fft.cpp
    fft-dipolar.cpp
    forcecap.cpp
    forces.cpp
    galilei.cpp
    ghmc.cpp
    ghosts.cpp
    global.cpp
    grid.cpp
    halo.cpp
    immersed_boundaries.cpp
    initialize.cpp
    integrate.cpp
    lattice.cpp
    layered.cpp
    metadynamics.cpp
    minimize_energy.cpp
    MpiCallbacks.cpp
    nemd.cpp
    npt.cpp
    nsquare.cpp
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
install(TARGETS EspressoConfig LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(EspressoConfig PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(EspressoConfig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(EspressoCore SHARED ${EspressoCore_SRC} config-version.cpp)
add_dependencies(EspressoCore EspressoConfig)
install(TARGETS EspressoCore LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
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

target_link_libraries(EspressoCore PRIVATE EspressoConfig ${LIBRARIES} Actor ObjectInFluid ImmersedBoundary Shapes Constraints Accumulators cluster_analysis)


if(SCAFACOS)
  target_link_libraries(EspressoCore PRIVATE Scafacos)
endif(SCAFACOS)

if(GSL_FOUND)
  target_include_directories(EspressoCore PUBLIC ${GSL_INCLUDE_DIRS})
  target_link_libraries(EspressoCore PRIVATE GSL::gsl GSL::gslcblas)
endif(GSL_FOUND)

if(CUDA)
    set(EspressoCuda_SRC
        cuda_common_cuda.cu
        cuda_init_cuda.cu
        EspressoSystemInterface_cuda.cu
    )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC})
  install(TARGETS EspressoCuda DESTINATION ${PYTHON_INSTDIR}/espressomd)
  set_target_properties(EspressoCuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(EspressoCuda EspressoConfig)

  target_link_libraries(EspressoCuda ${CUDA_CUFFT_LIBRARIES})
  target_link_libraries(EspressoCore PRIVATE EspressoCuda)
endif(CUDA)

# Subdirectories
add_subdirectory(observables)
target_link_libraries(EspressoCore PUBLIC Observables)
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
add_subdirectory(grid_based_algorithms)
target_link_libraries(EspressoCore PRIVATE lbboundaries_interface)
target_link_libraries(EspressoCore PRIVATE lbcpu)
if(CUDA)
  target_link_libraries(EspressoCore PRIVATE electrokinetics_core)
  target_link_libraries(EspressoCore PRIVATE fd_electrostatics)
  target_link_libraries(EspressoCore PRIVATE lbgpu)
endif(CUDA)
add_subdirectory(immersed_boundary)
add_subdirectory(virtual_sites)
target_link_libraries(EspressoCore PRIVATE VirtualSites)
if(CUDA)
  target_link_libraries(EspressoCore PRIVATE VirtualSitesCuda)
endif(CUDA)
add_subdirectory(nonbonded_interactions)
add_subdirectory(object-in-fluid)
add_subdirectory(shapes)
add_subdirectory(constraints)
add_subdirectory(cluster_analysis)
add_subdirectory(electrostatics_magnetostatics)
target_link_libraries(EspressoCore PRIVATE electrostatics_magnetostatics_core)
if(CUDA)
  target_link_libraries(EspressoCore PRIVATE electrostatics_magnetostatics_core_cuda)
endif()
