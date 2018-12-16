add_library(sd_interface SHARED sd_interface.cpp)
target_include_directories(sd_interface PRIVATE
  ${CMAKE_BINARY_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/core/stokesian_dynamics
  ${CMAKE_SOURCE_DIR}/libs/stokesian_dynamics
)
target_link_libraries(sd_interface PUBLIC EspressoConfig Boost::mpi MPI::MPI_CXX Random123 utils)

if (CUDA)
  target_link_libraries(sd_interface PUBLIC sd_gpu)
endif()

if(BLAS AND LAPACK)
  target_link_libraries(sd_interface PUBLIC sd_cpu)
endif()
