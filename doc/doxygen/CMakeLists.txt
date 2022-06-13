find_package(Doxygen)
if(DOXYGEN_FOUND)
  add_custom_command(
    OUTPUT doxy-features
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_doxyconfig.py
            ${CMAKE_SOURCE_DIR}/src/config doxy-features
    DEPENDS Espresso::config)

  set(DOXYGEN_BIB_IN ${CMAKE_SOURCE_DIR}/doc/bibliography.bib)
  set(DOXYGEN_BIB_OUT ${CMAKE_CURRENT_BINARY_DIR}/bibliography.bib)

  # transform BibTeX DOI fields into URL fields (bibliographic styles available
  # to Doxygen do not process the DOI field)
  add_custom_command(
    OUTPUT ${DOXYGEN_BIB_OUT}
    COMMAND
      sed -r
      "'s_^ *doi *= *([^0-9]+)(10\\.[0-9]+)_url=\\1https://dx.doi.org/\\2_'"
      ${DOXYGEN_BIB_IN} > ${DOXYGEN_BIB_OUT}
    DEPENDS ${DOXYGEN_BIB_IN})

  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

  add_custom_target(
    doxygen COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxy-features ${DOXYGEN_BIB_OUT}
    VERBATIM)

endif(DOXYGEN_FOUND)
