add_custom_target(dg
        COMMAND ${CMAKE_BINARY_DIR}/doc/latexit.sh
        ${CMAKE_CURRENT_SOURCE_DIR} dg
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${CMAKE_BINARY_DIR}/doc/latexit.sh
        ${CMAKE_SOURCE_DIR}/doc/logo/logo.pdf
                ${CMAKE_CURRENT_SOURCE_DIR}/todonotes.sty
)
