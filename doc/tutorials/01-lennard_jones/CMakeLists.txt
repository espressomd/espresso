configure_file(figures/lennard-jones-potential.png ${CMAKE_CURRENT_BINARY_DIR}/figures/lennard-jones-potential.png COPYONLY)
configure_file(01-lennard_jones.ipynb ${CMAKE_CURRENT_BINARY_DIR}/01-lennard_jones.ipynb COPYONLY)

add_custom_target(tutorials_01)

nb_export(TUTORIAL tutorials_01 FILE "01-lennard_jones.ipynb" HTML_RUN VAR_SUBST "sampling_iterations=4000")
