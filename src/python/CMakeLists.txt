


set(PYTHON_DIR ${CMAKE_CURRENT_BINARY_DIR})
configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso)

set(PYTHON_DIR ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR})
configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pypresso)
option(INSTALL_PYPRESSO "Install pypresso script, not needed when Espresso is installed in /usr" ON)
if(INSTALL_PYPRESSO)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pypresso
        DESTINATION ${BINDIR})
endif(INSTALL_PYPRESSO)

# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)

