configure_file(12-constant_pH.ipynb ${CMAKE_CURRENT_BINARY_DIR}/12-constant_pH.ipynb COPYONLY)

add_custom_target(tutorial_12)

nb_export(TUTORIAL tutorial_12 FILE "12-constant_pH.ipynb" HTML_RUN)

