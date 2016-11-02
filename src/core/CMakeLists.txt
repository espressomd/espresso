file(GLOB EspressoCore_SRC
          "*.cpp"
          )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-features.hpp ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  ${CMAKE_SOURCE_DIR}/src/features.def
  ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
  ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  )

add_custom_command(OUTPUT config-version.cpp
                   COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

# Custom target for everything that needs to be triggered when the myconfig changes
add_custom_target(EspressoConfig DEPENDS config-features.hpp config-features.cpp)
add_dependencies(EspressoConfig myconfig)

add_library(EspressoCore SHARED ${EspressoCore_SRC} config-features.cpp config-version.cpp)
set_target_properties(EspressoCore PROPERTIES SOVERSION ${SOVERSION})
install(TARGETS EspressoCore LIBRARY DESTINATION ${LIBDIR})
add_dependencies(EspressoCore EspressoConfig)

target_link_libraries(EspressoCore ${LIBRARIES} Actor ObjectInFluid ImmersedBoundary)

if(SCAFACOS)
  target_link_libraries(EspressoCore Scafacos)
endif(SCAFACOS)

set_target_properties(EspressoCore PROPERTIES MACOSX_RPATH TRUE)

if(CUDA)
  file(GLOB EspressoCuda_SRC
          "*.cu"
          )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC})
  install(TARGETS EspressoCuda DESTINATION ${LIBDIR})
  add_dependencies(EspressoCuda EspressoConfig)

  set_target_properties(EspressoCuda PROPERTIES MACOSX_RPATH TRUE)
  CUDA_ADD_CUFFT_TO_TARGET(EspressoCuda)
  target_link_libraries(EspressoCore EspressoCuda)
endif(CUDA)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)

add_subdirectory(actor)
add_subdirectory(immersed_boundary)
add_subdirectory(object-in-fluid)
if(SCAFACOS)
  add_subdirectory(scafacos)
endif(SCAFACOS)
