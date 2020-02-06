file(COPY "EXERCISES" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "SOLUTIONS" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "FIGURES" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(06-active_matter.ipynb ${CMAKE_CURRENT_BINARY_DIR}/06-active_matter.ipynb COPYONLY)

add_custom_target(tutorials_06)

nb_export(TUTORIAL tutorials_06 FILE "06-active_matter.ipynb" HTML_RUN)
