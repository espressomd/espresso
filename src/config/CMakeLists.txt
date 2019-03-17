configure_file(${CMAKE_SOURCE_DIR}/cmake/cmake_config.cmakein
               cmake_config.hpp
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/config-features.hpp ${CMAKE_CURRENT_BINARY_DIR}/config-features.cpp
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_featureconfig.py
  ${CMAKE_CURRENT_SOURCE_DIR}/features.def
  ${CMAKE_CURRENT_BINARY_DIR}/config-features.hpp
  ${CMAKE_CURRENT_BINARY_DIR}/config-features.cpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def ${CMAKE_CURRENT_SOURCE_DIR}/gen_featureconfig.py
  )
add_custom_target(
  generate_config_features
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/config-features.hpp ${CMAKE_CURRENT_BINARY_DIR}/config-features.cpp
)

add_custom_target(check_myconfig
  COMMAND ${PYTHON_EXECUTABLE}
  ${CMAKE_CURRENT_SOURCE_DIR}/check_myconfig.py
  ${CMAKE_CXX_COMPILER}
  ${CMAKE_CURRENT_SOURCE_DIR}/features.def
  ${CMAKE_CURRENT_BINARY_DIR}/myconfig-final.hpp
  cmake_config.hpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def
          myconfig
  )

execute_process(COMMAND ${PYTHON_EXECUTABLE}  
                ${CMAKE_CURRENT_SOURCE_DIR}/gen_sampleconfig.py 
                ${CMAKE_CURRENT_SOURCE_DIR}/features.def
                OUTPUT_FILE ${CMAKE_BINARY_DIR}/myconfig-sample.hpp
                )

add_library(EspressoConfig SHARED config-features.cpp)
add_dependencies(EspressoConfig myconfig check_myconfig generate_config_features)
install(TARGETS EspressoConfig LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(EspressoConfig PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(EspressoConfig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Parse repository info from git if available
include(version)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.hpp
  )
