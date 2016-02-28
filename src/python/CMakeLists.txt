configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso)
# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso DEPENDS ${CMAKE_BINARY_DIR}/pypresso)
# Make pypresso depend on the actual python module
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)
