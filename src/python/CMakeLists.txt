


set(PYTHON_DIR ${CMAKE_CURRENT_BINARY_DIR})
configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso)

set(PYTHON_DIR ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR})
configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pypresso)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pypresso
        DESTINATION ${BINDIR})

# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)

