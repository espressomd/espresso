if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -G")
else()
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3")
endif()
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30")

file(GLOB EspressoCore_SRC
          "*.cpp"
          )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT config-features.cpp config-features.hpp COMMAND
      ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
      ${CMAKE_SOURCE_DIR}/src/features.def
      ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
      ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
      DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_BINARY_DIR}/src/core/myconfig-final.hpp
)

add_custom_target(EspressoConfig ALL DEPENDS config-features.cpp config-features.hpp)

add_custom_command(OUTPUT config-version.cpp
                   COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

add_library(EspressoCore SHARED ${EspressoCore_SRC}
                        ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
                        ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)
target_link_libraries(EspressoCore ${LIBRARIES} Actor ObjectInFluid Electrostatics ImmersedBoundary)
set_target_properties(EspressoCore PROPERTIES MACOSX_RPATH TRUE)

if(WITH_CUDA)
  file(GLOB EspressoCuda_SRC
          "*.cu"
          )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC} config-version.cpp config-features.cpp)
  
  add_dependencies(EspressoCore EspressoCuda)
  set_target_properties(EspressoCuda PROPERTIES MACOSX_RPATH TRUE)
  CUDA_ADD_CUFFT_TO_TARGET(EspressoCuda)
  target_link_libraries(EspressoCore EspressoCuda)
endif(WITH_CUDA)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)

add_subdirectory(actor)
add_subdirectory(immersed_boundary)
add_subdirectory(object-in-fluid)
add_subdirectory(electrostatics)
