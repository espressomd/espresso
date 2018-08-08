set(EspressoCore_SRC
    accumulators.cpp
    bmhtf-nacl.cpp
    buckingham.cpp
    cells.cpp
    collision.cpp
    comfixed_global.cpp
    communication.cpp
    constraints.cpp
    cos2.cpp
    cuda_init.cpp
    cuda_interface.cpp
    debug.cpp
    debye_hueckel.cpp
    dihedral.cpp
    domain_decomposition.cpp
    dpd.cpp
    elc.cpp
    electrokinetics.cpp
    electrokinetics_pdb_parse.cpp
    energy.cpp
    errorhandling.cpp
    EspressoSystemInterface.cpp
    fene.cpp
    fft-common.cpp
    fft.cpp
    fft-dipolar.cpp
    forcecap.cpp
    forces.cpp
    galilei.cpp
    gaussian.cpp
    gb.cpp
    Geometry.cpp
    ghmc.cpp
    ghosts.cpp
    global.cpp
    grid.cpp
    halo.cpp
    harmonic_dumbbell.cpp
    hat.cpp
    hertzian.cpp
    hydrogen_bond.cpp
    iccp3m.cpp
    immersed_boundaries.cpp
    initialize.cpp
    integrate.cpp
    interaction_data.cpp
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
    morse.cpp
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
    quartic.cpp
    random.cpp
    rattle.cpp
    reaction_ensemble.cpp
    reaction_field.cpp
    rotate_system.cpp
    rotation.cpp
    RuntimeErrorCollector.cpp
    RuntimeError.cpp
    RuntimeErrorStream.cpp
    scafacos.cpp
    soft_sphere.cpp
    specfunc.cpp
    statistics_chain.cpp
    statistics_cluster.cpp
    statistics.cpp
    statistics_fluid.cpp
    steppot.cpp
    swimmer_reaction.cpp
    SystemInterface.cpp
    tab.cpp
    thermostat.cpp
    topology.cpp
    tuning.cpp
    twist_stack.cpp
    umbrella.cpp
    utils.cpp
    Vector.cpp
    virtual_sites.cpp
)


if( WITH_COVERAGE )
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -Og")
  if( CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
    link_libraries(gcov)
  endif()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-features.hpp ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  ${CMAKE_SOURCE_DIR}/src/features.def
  ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
  ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  )
add_custom_target(generate_config_features DEPENDS ${CMAKE_BINARY_DIR}/src/core/config-features.hpp)

add_custom_command(OUTPUT config-version.cpp
                   COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

# Custom target for everything that needs to be triggered when the myconfig changes
add_library(EspressoConfig SHARED config-features.cpp)
add_dependencies(EspressoConfig myconfig check_myconfig)
install(TARGETS EspressoConfig LIBRARY DESTINATION ${LIBDIR})
set_target_properties(EspressoConfig PROPERTIES SOVERSION ${SOVERSION})

add_library(EspressoCore SHARED ${EspressoCore_SRC} config-version.cpp)
add_dependencies(EspressoCore generate_config_features)
install(TARGETS EspressoCore LIBRARY DESTINATION ${LIBDIR})
set_target_properties(EspressoCore PROPERTIES SOVERSION ${SOVERSION})

target_link_libraries(EspressoCore EspressoConfig ${LIBRARIES} Actor ObjectInFluid ImmersedBoundary Shapes Constraints Observables Accumulators cluster_analysis VirtualSites)

if(SCAFACOS)
  target_link_libraries(EspressoCore Scafacos)
endif(SCAFACOS)
if(GSL)
  target_include_directories(EspressoCore PUBLIC ${GSL_INCLUDE_DIRS})
endif(GSL)


if(GSL)
  target_link_libraries(EspressoCore gsl gslcblas)
endif(GSL)

if(GSL)
  target_link_libraries(EspressoCore gsl gslcblas)
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
  target_link_libraries(EspressoCore EspressoCuda)
endif(CUDA)

# Subdirectories
add_subdirectory(observables)
add_subdirectory(accumulators)
add_subdirectory(io)
if(HDF5_FOUND)
    target_link_libraries(EspressoCore H5mdCore)
endif(HDF5_FOUND)
target_link_libraries(EspressoCore mpiio)
target_link_libraries(EspressoCore pdbreader)
target_link_libraries(EspressoCore bonded_interactions)
target_link_libraries(EspressoCore nonbonded_interactions)

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
