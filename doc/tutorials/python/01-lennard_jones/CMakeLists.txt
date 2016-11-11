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


# FIXME here: for the          XX  give it a proper number (01, 02...)
add_custom_target(tutorials_py_01 DEPENDS ${BASENAME}.pdf)

configure_file(scripts/lj_tutorial.py ${CMAKE_CURRENT_BINARY_DIR}/lj_tutorial.py COPYONLY)
configure_file(scripts/msd.py ${CMAKE_CURRENT_BINARY_DIR}/msd.py COPYONLY)
configure_file(scripts/two-component.py ${CMAKE_CURRENT_BINARY_DIR}/two-component.py COPYONLY)
configure_file(scripts/two-component-visualization.py ${CMAKE_CURRENT_BINARY_DIR}/two-component-visualization.py COPYONLY)
