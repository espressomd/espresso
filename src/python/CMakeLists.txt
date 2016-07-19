configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso)
option(INSTALL_PYPRESSO "Install pypresso files, no need when installing under /usr" ON)
if(INSTALL_PYPRESSO)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/pypresso DESTINATION bin)
endif(INSTALL_PYPRESSO)
# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)
