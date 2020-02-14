configure_file(08-visualization.ipynb ${CMAKE_CURRENT_BINARY_DIR}/08-visualization.ipynb COPYONLY)

add_custom_target(tutorial_08)

nb_export(TUTORIAL tutorial_08 FILE "08-visualization.ipynb" HTML_RUN)
