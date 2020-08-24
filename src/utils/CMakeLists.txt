add_library(EspressoUtils INTERFACE)
target_include_directories(
  EspressoUtils INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                          $<INSTALL_INTERFACE:include>)
target_link_libraries(EspressoUtils INTERFACE Boost::serialization Boost::mpi
                                              MPI::MPI_CXX)

install(TARGETS EspressoUtils LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
