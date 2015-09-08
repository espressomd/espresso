
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -gencode arch=compute_35,code=sm_35)

file(GLOB_RECURSE EspressoCore_SRC
          "*.cpp"
          )

file(GLOB_RECURSE EspressoCuda_SRC
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

cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC})

add_dependencies(EspressoCuda
                 config-features
)

add_library(EspressoCore SHARED ${EspressoCore_SRC}
                        ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
                        ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

add_dependencies(EspressoCore
                 EspressoCuda
                 config-features
                 )

target_link_libraries(EspressoCore EspressoCuda)

add_subdirectory(actor)
