configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso)
configure_file(pypresso_install.cmake.cmakein pypresso_install.cmake)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/pypresso_install.cmake)
# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)

