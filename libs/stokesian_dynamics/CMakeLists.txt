if(STOKESIAN_DYNAMICS_GPU)
  cuda_add_library(sd_gpu SHARED sd_gpu.cu)
  target_link_libraries(sd_gpu PRIVATE
    ${CUDA_CUBLAS_LIBRARIES}
    ${CUDA_cusolver_LIBRARY})
  target_compile_definitions(sd_gpu PUBLIC USE_THRUST)
  install(TARGETS sd_gpu LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
endif()

if(STOKESIAN_DYNAMICS)
  add_library(sd_cpu SHARED sd_cpu.cpp)
  target_link_libraries(sd_cpu 
    PRIVATE
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      Boost::boost
      Random123)
  if(STOKESIAN_DYNAMICS_GPU)
    target_compile_definitions(sd_cpu PUBLIC USE_THRUST)
  endif()
  install(TARGETS sd_cpu LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
endif()
