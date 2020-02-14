configure_file(12-constant_pH.ipynb ${CMAKE_CURRENT_BINARY_DIR}/12-constant_pH.ipynb COPYONLY)

add_custom_target(tutorial_12)
add_dependencies(tutorials tutorial_12)

nb_export(TARGET tutorial_12 SUFFIX "" FILE "12-constant_pH.ipynb" HTML_RUN)

