set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3 -arch=sm_20")

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


add_custom_command(OUTPUT config-features.cpp config-features.hpp COMMAND
      ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
      ${CMAKE_SOURCE_DIR}/src/features.def
      ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
      ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
      DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_BINARY_DIR}/src/core/myconfig-final.hpp
)

add_custom_command(OUTPUT config-version.cpp
                   COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC} config-version.cpp config-features.cpp)


add_library(EspressoCore SHARED ${EspressoCore_SRC}
                        ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
                        ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

add_dependencies(EspressoCore
                 EspressoCuda
                 )

target_link_libraries(EspressoCore EspressoCuda)

add_subdirectory(actor)
