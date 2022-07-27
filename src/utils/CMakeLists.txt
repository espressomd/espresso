add_library(Espresso_utils INTERFACE)
add_library(Espresso::utils ALIAS Espresso_utils)
add_library(Espresso_utils_mpi INTERFACE)
add_library(Espresso::utils::mpi ALIAS Espresso_utils)
target_include_directories(
  Espresso_utils
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(Espresso_utils INTERFACE Boost::serialization)
target_link_libraries(Espresso_utils_mpi INTERFACE Espresso::utils Boost::mpi
                                                   MPI::MPI_CXX)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
