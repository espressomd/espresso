# temporarily disabled to at least be able to test the CPU version
# if(STOKESIAN_DYNAMICS_GPU)
#   cuda_add_library(sd_gpu SHARED sd_gpu.cu)
#   target_link_libraries(sd_gpu PRIVATE
#     ${CUDA_CUBLAS_LIBRARIES}
#     ${CUDA_cusolver_LIBRARY})
# endif()

if(STOKESIAN_DYNAMICS)
  add_library(sd_cpu SHARED sd_cpu.cpp)
  target_link_libraries(sd_cpu PRIVATE
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES})
endif()
