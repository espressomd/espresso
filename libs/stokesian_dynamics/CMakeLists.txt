if(STOKESIAN_DYNAMICS_GPU)
  add_gpu_library(sd_gpu SHARED sd_gpu.cu)
  target_link_libraries(sd_gpu PRIVATE
    ${CUDA_CUBLAS_LIBRARIES}
    ${CUDA_cusolver_LIBRARY})
  target_compile_definitions(sd_gpu PRIVATE SD_USE_THRUST SD_USE_GPU)
  if(HIP)
    target_compile_definitions(sd_gpu PRIVATE SD_USE_ROCM)
  else()
    target_compile_definitions(sd_gpu PRIVATE SD_USE_CUDA)
  endif()
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

  # In case the GPU is used, Thrust is present and can be used to parallelize
  # the CPU code, too. The standard compiler needs to be told the location of
  # Thrust
  if(STOKESIAN_DYNAMICS_GPU)
    target_compile_definitions(sd_cpu PUBLIC SD_USE_THRUST PRIVATE THRUST_DEVICE_SYSTEM=4)
    #target_compile_definitions(sd_cpu PUBLIC USE_THRUST)
    target_include_directories(sd_cpu PRIVATE 
      ${CUDA_INCLUDE_DIRS}
      ${HIP_ROOT_DIR}/include
      ${ROCM_HOME}/include)
  endif()
  install(TARGETS sd_cpu LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
endif()
