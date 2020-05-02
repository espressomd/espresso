add_library(ScriptInterface SHARED initialize.cpp ScriptInterfaceBase.cpp
                                   ParallelScriptInterface.cpp)

add_subdirectory(accumulators)
add_subdirectory(collision_detection)
add_subdirectory(constraints)
add_subdirectory(cluster_analysis)
add_subdirectory(lbboundaries)
add_subdirectory(virtual_sites)
add_subdirectory(observables)
add_subdirectory(pair_criteria)
add_subdirectory(mpiio)
add_subdirectory(shapes)
add_subdirectory(h5md)

install(TARGETS ScriptInterface
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

target_link_libraries(
  ScriptInterface PRIVATE EspressoConfig EspressoCore
  PUBLIC mpiio utils MPI::MPI_CXX shapes core_cluster_analysis
  PRIVATE cxx_interface)

target_include_directories(ScriptInterface PUBLIC ${CMAKE_SOURCE_DIR}/src)
