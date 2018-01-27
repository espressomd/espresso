get_filename_component(BASENAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# we assume the tex filename is the same as the directory name
add_custom_command( 
  OUTPUT ${BASENAME}.pdf
  COMMAND sh ../../latexit.sh
          ${CMAKE_CURRENT_SOURCE_DIR}:${CMAKE_CURRENT_SOURCE_DIR}/../common
          ${BASENAME}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.tex
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/../common
)


# FIXME here: for the          XX  give it a proper number (01, 02...)
add_custom_target(tutorials_02 DEPENDS ${BASENAME}.pdf)

configure_file(scripts/nacl.py ${CMAKE_CURRENT_BINARY_DIR}/nacl.py COPYONLY)
configure_file(scripts/nacl_units_confined.py ${CMAKE_CURRENT_BINARY_DIR}/nacl_units_confined.py COPYONLY)
configure_file(scripts/nacl_units_confined_vis.py ${CMAKE_CURRENT_BINARY_DIR}/nacl_units_confined_vis.py COPYONLY)
configure_file(scripts/nacl_units.py ${CMAKE_CURRENT_BINARY_DIR}/nacl_units.py COPYONLY)
configure_file(scripts/nacl_units_vis.py ${CMAKE_CURRENT_BINARY_DIR}/nacl_units_vis.py COPYONLY)
