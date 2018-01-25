get_filename_component(BASENAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_custom_command( 
  OUTPUT ${BASENAME}.pdf
  COMMAND sh ../../../latexit.sh 
          ${CMAKE_CURRENT_SOURCE_DIR}:${CMAKE_CURRENT_SOURCE_DIR}/../../common
          ${BASENAME}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.tex
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/../../common
)

add_custom_target(tutorials_py_06 DEPENDS ${BASENAME}.pdf)
