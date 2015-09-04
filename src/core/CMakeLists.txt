
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -gencode arch=compute_35,code=sm_35)

file(GLOB_RECURSE libEspresso_SRC
          "*.cpp"
          )

file(GLOB_RECURSE libEspresso_cuda_SRC
          "*.cu"
          )

include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)

cuda_include_directories(${CMAKE_SOURCE_DIR}/src/core)
cuda_include_directories(${CMAKE_BINARY_DIR}/src/core)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-features.cpp ${CMAKE_BINARY_DIR}/src/core/config-features.hpp COMMAND
      ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
      ${CMAKE_SOURCE_DIR}/src/features.def
      ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
      ${CMAKE_BINARY_DIR}/src/core/config-features.cpp DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def
)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

cuda_add_library(libEspresso_cuda ${libEspresso_cuda_SRC} ${CMAKE_BINARY_DIR}/src/core/config-features.cpp)

add_library(libEspresso STATIC ${libEspresso_SRC}
${CMAKE_BINARY_DIR}/src/core/config-features.cpp
${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)
target_link_libraries(libEspresso libEspresso_cuda)

add_subdirectory(actor)
