configure_file(08-visualization.ipynb ${CMAKE_CURRENT_BINARY_DIR}/08-visualization.ipynb COPYONLY)

add_custom_target(tutorial_08)
add_dependencies(tutorials tutorial_08)

nb_export(TARGET tutorial_08 SUFFIX "" FILE "08-visualization.ipynb" HTML_RUN)
