if (CUDA)
  cuda_add_library(sd_gpu SHARED sd_gpu.cu)
  target_link_libraries(sd_gpu
    ${CUDA_CUBLAS_LIBRARIES}
    ${CUDA_cusolver_LIBRARY})
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_library(sd_cpu SHARED sd_cpu.cpp)
target_link_libraries(sd_cpu
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES})
