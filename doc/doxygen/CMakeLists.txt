find_package(Doxygen)
if(DOXYGEN_FOUND)
  add_custom_command(OUTPUT doxy-features
    COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/gen_doxyconfig.py
    ${CMAKE_SOURCE_DIR}/src/config doxy-features
    DEPENDS EspressoConfig
  )

  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

  add_custom_target(doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxy-features
    VERBATIM
  )

endif(DOXYGEN_FOUND)
