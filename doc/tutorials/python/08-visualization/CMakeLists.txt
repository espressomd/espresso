get_filename_component(BASENAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# we assume the tex filename is the same as the directory name
add_custom_command( 
  OUTPUT ${BASENAME}.pdf
  COMMAND sh ../../../latexit.sh 
          ${CMAKE_CURRENT_SOURCE_DIR}:${CMAKE_CURRENT_SOURCE_DIR}/../../common
          ${BASENAME}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.tex
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/../../common
)


add_custom_target(tutorials_py_08 DEPENDS ${BASENAME}.pdf)

configure_file(scripts/plotting.py ${CMAKE_CURRENT_BINARY_DIR}/plotting.py COPYONLY)
configure_file(scripts/simulation.py ${CMAKE_CURRENT_BINARY_DIR}/simulation.py COPYONLY)
configure_file(scripts/visualization.py ${CMAKE_CURRENT_BINARY_DIR}/visualization.py COPYONLY)
