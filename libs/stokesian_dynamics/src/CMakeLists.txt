# TODO: License header

if(STOKESIAN_DYNAMICS_GPU)
  add_gpu_library(sd_gpu SHARED sd_gpu.cu)
  add_library(StokesianDynamics::sd_gpu ALIAS sd_gpu)
  target_link_libraries(sd_gpu PRIVATE
    stokesian_dynamics
    ${CUDA_CUBLAS_LIBRARIES}
    ${CUDA_cusolver_LIBRARY})
  target_compile_definitions(sd_gpu PRIVATE SD_USE_THRUST SD_USE_GPU)
  if(HIP_VERSION)
    target_compile_definitions(sd_gpu PRIVATE SD_USE_ROCM)
  else()
    target_compile_definitions(sd_gpu PRIVATE SD_USE_CUDA)
  endif()
  install(TARGETS sd_gpu
    EXPORT stokesiandynamics-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  set_target_properties(sd_gpu PROPERTIES EXPORT_NAME sd_gpu)
endif()

if(STOKESIAN_DYNAMICS)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(Boost 1.65 REQUIRED)
  add_library(sd_cpu SHARED sd_cpu.cpp)
  add_library(StokesianDynamics::sd_cpu ALIAS sd_cpu)
  target_link_libraries(sd_cpu
    PRIVATE
      stokesian_dynamics
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      Boost::boost
      Random123)

  # In case the GPU is used, Thrust is present and can be used to parallelize
  # the CPU code, too. The standard compiler needs to be told the location of
  # Thrust
  if(STOKESIAN_DYNAMICS_GPU)
    target_compile_definitions(sd_cpu PUBLIC SD_USE_THRUST PRIVATE THRUST_DEVICE_SYSTEM=4)
    if(HIP_VERSION)
      target_include_directories(sd_cpu PRIVATE
        ${HIP_ROOT_DIR}/include
        ${ROCM_HOME}/include)
    else()
      target_include_directories(sd_cpu PRIVATE
        ${CUDA_INCLUDE_DIRS})
    endif()
  endif()
  install(TARGETS sd_cpu
    EXPORT stokesiandynamics-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  set_target_properties(sd_cpu PROPERTIES EXPORT_NAME sd_cpu)
endif()
