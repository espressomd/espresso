find_program( IPYTHON_EXECUTABLE NAMES ipython )

# Configure pypresso for build dir
set(PYTHON_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(PYTHON_FRONTEND ${PYTHON_EXECUTABLE})
configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso @ONLY)

if( IPYTHON_EXECUTABLE )
  set(PYTHON_FRONTEND ${IPYTHON_EXECUTABLE})
  configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/ipypresso @ONLY)
endif( )

# Configure pypresso for install dir
set(PYTHON_DIR ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR})

set(PYTHON_FRONTEND ${PYTHON_EXECUTABLE})
configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pypresso @ONLY)

if( IPYTHON_EXECUTABLE )
  set(PYTHON_FRONTEND ${IPYTHON_EXECUTABLE})
  configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/ipypresso @ONLY)
endif( )


option(INSTALL_PYPRESSO "Install pypresso script, not needed when Espresso is installed in /usr" ON)
if(INSTALL_PYPRESSO)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pypresso
        DESTINATION ${BINDIR})

  if( IPYTHON_EXECUTABLE )
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ipypresso
          DESTINATION ${BINDIR})
  endif( )
endif(INSTALL_PYPRESSO)

# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd)

add_subdirectory(espressomd)
