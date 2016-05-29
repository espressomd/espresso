configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso)
# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)
