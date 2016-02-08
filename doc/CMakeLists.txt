add_subdirectory(doxygen)
add_subdirectory(logo)
add_subdirectory(ug)
add_subdirectory(dg)

find_program(PDFLATEX pdflatex)
find_program(MAKEINDEX makeindex)

set(abs_srcdir ${CMAKE_CURRENT_SOURCE_DIR})
set(abs_top_srcdir ${CMAKE_SOURCE_DIR})

configure_file(latexit.sh.in
               ${CMAKE_BINARY_DIR}/doc/latexit.sh
               @ONLY
              )

add_custom_target(doc DEPENDS logo ug dg doxygen)
