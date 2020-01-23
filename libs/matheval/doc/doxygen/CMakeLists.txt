find_package(Doxygen)

if(DOXYGEN_FOUND)
  add_custom_target(doxygen
    COMMAND
    SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/../../
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
    COMMENT "Building HTML documentation with Doxygen")
endif()
