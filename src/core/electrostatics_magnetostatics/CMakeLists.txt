add_library(electrostatics_magnetostatics_core SHARED
    debye_hueckel.cpp
    elc.cpp
    icc.cpp
    maggs.cpp
    magnetic_non_p3m_methods.cpp
    mdlc_correction.cpp
    mmm1d.cpp
    mmm2d.cpp
    mmm-common.cpp
    p3m-common.cpp
    p3m.cpp
    p3m-dipolar.cpp
    p3m_gpu.cpp
    scafacos.cpp
)

if(CUDA)
    set(electrostatics_core_cuda_src
        p3m_gpu_cuda.cu
        p3m_gpu_error_cuda.cu
    )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  cuda_add_library(electrostatics_magnetostatics_core_cuda SHARED ${electrostatics_core_cuda_src})
  install(TARGETS electrostatics_magnetostatics_core_cuda DESTINATION ${LIBDIR})
  set_target_properties(electrostatics_magnetostatics_core_cuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(electrostatics_magnetostatics_core_cuda EspressoConfig)

  target_link_libraries(electrostatics_magnetostatics_core_cuda ${CUDA_CUFFT_LIBRARIES})
  target_link_libraries(electrostatics_magnetostatics_core_cuda EspressoCuda)
endif(CUDA)

target_include_directories(electrostatics_magnetostatics_core PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
target_compile_options(electrostatics_magnetostatics_core PRIVATE -fPIC)
add_dependencies(electrostatics_magnetostatics_core EspressoConfig)

if(SCAFACOS)
  add_subdirectory(scafacos)
endif(SCAFACOS)
