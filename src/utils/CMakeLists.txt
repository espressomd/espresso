add_library(espresso_utils INTERFACE)
add_library(espresso::utils ALIAS espresso_utils)
add_library(espresso_utils_mpi INTERFACE)
add_library(espresso::utils::mpi ALIAS espresso_utils)
target_include_directories(
  espresso_utils
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(espresso_utils INTERFACE Boost::serialization)
target_link_libraries(espresso_utils_mpi INTERFACE espresso::utils Boost::mpi
                                                   MPI::MPI_CXX)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
