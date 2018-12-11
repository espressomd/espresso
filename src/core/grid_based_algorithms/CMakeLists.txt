add_subdirectory(lbboundaries)

if (CUDA)
  cuda_add_library(electrokinetics_core_cuda SHARED
    electrokinetics_cuda.cu
  )
  target_link_libraries(electrokinetics_core_cuda EspressoConfig ${CUDA_CUFFT_LIBRARIES})


  add_library(electrokinetics_core SHARED
    electrokinetics.cpp
    electrokinetics_pdb_parse.cpp
  )
  target_include_directories(electrokinetics_core PRIVATE
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_BINARY_DIR}/src/core
    )

  target_link_libraries(electrokinetics_core PRIVATE EspressoConfig electrokinetics_core_cuda)

  cuda_add_library(lbgpu_cuda SHARED
    lbgpu_cuda.cu
  )
  target_link_libraries(lbgpu_cuda EspressoConfig)

  add_library(lbgpu SHARED
    lbgpu.cpp
    lb-d3q19.cpp
  )
  target_include_directories(lbgpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_BINARY_DIR}/src/core
  )
  target_link_libraries(lbgpu PRIVATE lbgpu_cuda EspressoConfig)

  cuda_add_library(fd_electrostatics SHARED
    fd-electrostatics_cuda.cu
  )
  target_link_libraries(fd_electrostatics EspressoConfig ${CUDA_CUFFT_LIBRARIES})
endif()
