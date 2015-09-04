file(GLOB ug_TEXFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.tex)
file(GLOB ug_BIBFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.bib)
file(GLOB ug_FIGS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} figures/*.pdf)
file(GLOB ug_TEXFILES_abs *tex)

find_program(AWK awk)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/quickref.inp
                   COMMAND ${AWK} -f
                   ${CMAKE_CURRENT_SOURCE_DIR}/assemble_quickref.awk
                   ${ug_TEXFILES_abs} > ${CMAKE_CURRENT_BINARY_DIR}/quickref.inp
)

add_custom_target(ug
        COMMAND ${CMAKE_BINARY_DIR}/doc/latexit.sh
        ${CMAKE_CURRENT_SOURCE_DIR} ug
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/quickref.inp
                ${CMAKE_BINARY_DIR}/doc/latexit.sh
                ${CMAKE_BINARY_DIR}/doc/logo/logo.pdf
)

#add_latex_document(ug.tex
#        INPUTS ${ug_TEXFILES}
#        IMAGES ${ug_FIGS}
#        BIBFILES citations.bib otherrefs.bib
#        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/quickref.inp
#        USE_INDEX
#)

