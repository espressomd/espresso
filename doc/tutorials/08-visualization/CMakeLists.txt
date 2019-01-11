configure_file(08-visualization.ipynb ${CMAKE_CURRENT_BINARY_DIR}/08-visualization.ipynb COPYONLY)

add_custom_target(tutorials_08)

html_export(TUTORIAL tutorials_08 FILE "08-visualization.ipynb")
python_export(TUTORIAL tutorials_08 FILE "08-visualization.ipynb")
