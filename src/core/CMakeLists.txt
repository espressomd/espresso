
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -gencode arch=compute_35,code=sm_35)

file(GLOB_RECURSE libEspresso_SRC
          "*.cpp"
          )

file(GLOB_RECURSE libEspresso_cuda_SRC
          "*.cu"
          )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(config-features)
add_dependencies(config-features myconfig)

add_custom_command(TARGET config-features COMMAND
      ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
      ${CMAKE_SOURCE_DIR}/src/features.def
      ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
      ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
      DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def
)

add_custom_command(TARGET config-features
                   COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

cuda_add_library(libEspresso_cuda ${libEspresso_cuda_SRC})

add_dependencies(libEspresso_cuda 
                 config-features
)

add_library(libEspresso STATIC ${libEspresso_SRC}
                        ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
                        ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

add_dependencies(libEspresso
                 libEspresso_cuda
                 config-features
                 )

target_link_libraries(libEspresso libEspresso_cuda)

add_subdirectory(actor)
