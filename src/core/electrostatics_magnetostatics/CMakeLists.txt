if(CUDA)
    set(electrostatics_core_cuda_src
        p3m_gpu_cuda.cu
        p3m_gpu_error_cuda.cu
    )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})
  cuda_include_directories(${CMAKE_SOURCE_DIR}/src/config)
  cuda_include_directories(${CMAKE_BINARY_DIR}/src/config)

  cuda_add_library(electrostatics_magnetostatics_core_cuda SHARED ${electrostatics_core_cuda_src})
  install(TARGETS electrostatics_magnetostatics_core_cuda DESTINATION ${PYTHON_INSTDIR}/espressomd)
  set_target_properties(electrostatics_magnetostatics_core_cuda PROPERTIES SOVERSION ${SOVERSION})
  target_link_libraries(electrostatics_magnetostatics_core_cuda EspressoConfig)

  target_link_libraries(electrostatics_magnetostatics_core_cuda ${CUDA_CUFFT_LIBRARIES})
  target_link_libraries(electrostatics_magnetostatics_core_cuda EspressoCuda)
endif(CUDA)

if(SCAFACOS)
  add_subdirectory(scafacos)
endif(SCAFACOS)
