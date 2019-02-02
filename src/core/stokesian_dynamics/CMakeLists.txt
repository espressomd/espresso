add_library(sd_interface SHARED sd_interface.cpp)
target_include_directories(sd_interface PRIVATE
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_BINARY_DIR}/src/core
)
target_link_libraries(sd_interface PUBLIC EspressoConfig Boost::mpi MPI::MPI_CXX)

if (CUDA)
  cuda_add_library(sd_gpu SHARED sd_gpu.cu)
  target_link_libraries(sd_gpu PRIVATE
    ${CUDA_CUBLAS_LIBRARIES}
    ${CUDA_cusolver_LIBRARY})
  target_link_libraries(sd_interface PUBLIC sd_gpu)
endif()

if(BLAS AND LAPACK)
  add_library(sd_cpu SHARED sd_cpu.cpp)
  target_link_libraries(sd_cpu PRIVATE
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES})
  target_link_libraries(sd_interface PUBLIC sd_cpu)
endif()
