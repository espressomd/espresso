# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_SOURCE_DIR}/pypresso.in
  COMMAND ${CMAKE_COMMAND} -E configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}}/pypresso) )
# Make pypresso depend on the actual python module
add_dependencies(pypresso _init)

add_subdirectory(espressomd)
