add_library(
  Espresso_script_interface SHARED
  initialize.cpp ObjectHandle.cpp object_container_mpi_guard.cpp
  GlobalContext.cpp ContextManager.cpp ParallelExceptionHandler.cpp)
add_library(Espresso::script_interface ALIAS Espresso_script_interface)

add_subdirectory(accumulators)
add_subdirectory(bond_breakage)
add_subdirectory(cell_system)
add_subdirectory(collision_detection)
add_subdirectory(constraints)
add_subdirectory(cluster_analysis)
add_subdirectory(electrostatics)
add_subdirectory(interactions)
add_subdirectory(lbboundaries)
add_subdirectory(lees_edwards)
add_subdirectory(magnetostatics)
add_subdirectory(virtual_sites)
add_subdirectory(observables)
add_subdirectory(pair_criteria)
add_subdirectory(mpiio)
add_subdirectory(scafacos)
add_subdirectory(shapes)
add_subdirectory(h5md)
add_subdirectory(reaction_methods)

install(TARGETS Espresso_script_interface
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

target_link_libraries(
  Espresso_script_interface PRIVATE Espresso::config Espresso::core
  PUBLIC Espresso::utils MPI::MPI_CXX Espresso::shapes
  PRIVATE Espresso::cpp_flags)

target_include_directories(Espresso_script_interface
                           PUBLIC ${CMAKE_SOURCE_DIR}/src)

add_subdirectory(tests)
