add_library(Espresso_utils INTERFACE)
add_library(Espresso::utils ALIAS Espresso_utils)
target_include_directories(
  Espresso_utils
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)
target_link_libraries(Espresso_utils INTERFACE Boost::serialization Boost::mpi
                                               MPI::MPI_CXX)

install(TARGETS Espresso_utils LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
