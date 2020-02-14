file(COPY "EXERCISES" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "SOLUTIONS" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "FIGURES" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(06-active_matter.ipynb ${CMAKE_CURRENT_BINARY_DIR}/06-active_matter.ipynb COPYONLY)

add_custom_target(tutorial_06)
add_dependencies(tutorials tutorial_06)

nb_export(TARGET tutorial_06 SUFFIX "" FILE "06-active_matter.ipynb" HTML_RUN)
