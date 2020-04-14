find_package(Repa)

set(GenericDD_SRC generic_dd.cpp metric.cpp)

add_library(GenericDD SHARED ${GenericDD_SRC})
install(TARGETS GenericDD LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
add_dependencies(GenericDD EspressoConfig)

target_link_libraries(
  GenericDD PRIVATE EspressoConfig Profiler utils Boost::serialization
                    Boost::mpi MPI::MPI_CXX)

target_include_directories(GenericDD PRIVATE ${CMAKE_SOURCE_DIR}/src/core
                                             ${CMAKE_BINARY_DIR}/src/core)

if(REPA_FOUND)
  target_link_libraries(GenericDD PUBLIC ${REPA_LIBRARIES})
  target_include_directories(GenericDD PRIVATE ${REPA_INCLUDE_DIR})
  target_compile_definitions(GenericDD PRIVATE HAVE_REPA)
endif(REPA_FOUND)
