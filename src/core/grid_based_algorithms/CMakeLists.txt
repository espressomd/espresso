add_subdirectory(lbboundaries)

if (CUDA)
  cuda_add_library(electrokinetics_core SHARED
    electrokinetics.cpp
    electrokinetics_cuda.cu
    electrokinetics_pdb_parse.cpp
  )
  target_include_directories(electrokinetics_core PRIVATE
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_BINARY_DIR}/src/core
  )
  add_dependencies(electrokinetics_core EspressoConfig)
  target_link_libraries(electrokinetics_core "cufft")

  cuda_add_library(lbgpu SHARED
    lbgpu.cpp
    lb-d3q19.cpp
    lbgpu_cuda.cu
  )
  target_include_directories(lbgpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_BINARY_DIR}/src/core
  )
  add_dependencies(lbgpu EspressoConfig)

  cuda_add_library(fd_electrostatics SHARED
    fd-electrostatics_cuda.cu
  )
  add_dependencies(fd_electrostatics EspressoConfig)
  target_link_libraries(fd_electrostatics "cufft")
endif()

add_library(lbboundaries_interface SHARED
  lbboundaries.cpp
)
target_include_directories(lbboundaries_interface PRIVATE
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_BINARY_DIR}/src/core
)
add_dependencies(lbboundaries_interface EspressoConfig)

add_library(lbcpu SHARED
  lb.cpp
  lb-d3q19.cpp
)
target_include_directories(lbcpu PRIVATE
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_BINARY_DIR}/src/core
)
add_dependencies(lbcpu EspressoConfig)
