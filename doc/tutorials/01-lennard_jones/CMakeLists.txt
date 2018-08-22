configure_file(figures/lennard-jones-potential.png ${CMAKE_CURRENT_BINARY_DIR}/figures/lennard-jones-potential.png COPYONLY)
configure_file(01-lennard_jones.ipynb ${CMAKE_CURRENT_BINARY_DIR}/01-lennard_jones.ipynb COPYONLY)

add_custom_target(tutorials_01)

html_export(TUTORIAL tutorials_01 FILE "01-lennard_jones.ipynb")
python_export(TUTORIAL tutorials_01 FILE "01-lennard_jones.ipynb")
