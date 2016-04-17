set(logo_FILES
	logo.svg
	logo.pdf logo_500x500.png
	logo-small.svg
	logo-small.pdf logo_32x32.png logo_48x48.png logo_100x100.png
	logo-animated-100.gif logo-animated-200.gif logo-animated-500.gif
)


foreach(lf in logo_FILES)
        add_custom_command(OUTPUT ${lf}
                           COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${lf} ${CMAKE_CURRENT_BINARY_DIR}
                          )
endforeach()



