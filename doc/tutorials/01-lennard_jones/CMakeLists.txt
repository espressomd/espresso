configure_file(figures/lennard-jones-potential.png ${CMAKE_CURRENT_BINARY_DIR}/figures/lennard-jones-potential.png COPYONLY)
configure_file(01-lennard_jones.ipynb ${CMAKE_CURRENT_BINARY_DIR}/01-lennard_jones.ipynb COPYONLY)

add_custom_target(tutorial_01)
add_dependencies(tutorials tutorial_01)

nb_export(TARGET tutorial_01 SUFFIX "" FILE "01-lennard_jones.ipynb" HTML_RUN VAR_SUBST)
