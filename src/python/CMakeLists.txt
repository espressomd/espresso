# Configure pypresso for build dir
set(PYTHON_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(PYTHON_FRONTEND ${PYTHON_EXECUTABLE})
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(
    COMMAND "${PYTHON_FRONTEND}" "-c"
            "import sysconfig; print(sysconfig.get_config_var('exec_prefix'))"
    OUTPUT_VARIABLE DARWIN_EXEC_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PYTHON_FRONTEND_TMP
      "${DARWIN_EXEC_PREFIX}/Resources/Python.app/Contents/MacOS/Python")
  if(EXISTS ${PYTHON_FRONTEND_TMP})
    set(PYTHON_FRONTEND ${PYTHON_FRONTEND_TMP})
  endif()
endif()
configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso @ONLY)

if(IPYTHON_EXECUTABLE)
  set(PYTHON_FRONTEND ${IPYTHON_EXECUTABLE})
  configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/ipypresso @ONLY)
endif()

# Configure pypresso for install dir
set(PYTHON_DIR ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR})

set(PYTHON_FRONTEND ${PYTHON_EXECUTABLE})
configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pypresso @ONLY)

if(IPYTHON_EXECUTABLE)
  set(PYTHON_FRONTEND ${IPYTHON_EXECUTABLE})
  configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/ipypresso @ONLY)
endif()

option(INSTALL_PYPRESSO
       "Install pypresso script, not needed when ESPResSo is installed in /usr"
       ON)
if(INSTALL_PYPRESSO)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pypresso
          DESTINATION ${CMAKE_INSTALL_BINDIR})

  if(IPYTHON_EXECUTABLE)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ipypresso
            DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif(INSTALL_PYPRESSO)

# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd object_in_fluid)

add_subdirectory(espressomd)
add_subdirectory(object_in_fluid)
