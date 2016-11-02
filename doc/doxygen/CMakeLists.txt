find_package(Doxygen)
if(DOXYGEN_FOUND)

add_custom_command(OUTPUT doxy-features
                   COMMAND ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/gen_doxyconfig.py
                   ${CMAKE_SOURCE_DIR}/src/features.def doxy-features
                   DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def
)

add_custom_target(doxygen
  COMMAND BUILDDIR=${CMAKE_CURRENT_BINARY_DIR} SRCDIR=${CMAKE_SOURCE_DIR}
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxy-features
  COMMENT "Generating API documentation with Doxygen" VERBATIM
)

endif(DOXYGEN_FOUND)
