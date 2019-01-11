configure_file(scripts/plotting.py ${CMAKE_CURRENT_BINARY_DIR}/plotting.py COPYONLY)
configure_file(scripts/simulation.py ${CMAKE_CURRENT_BINARY_DIR}/simulation.py COPYONLY)
configure_file(scripts/visualization.py ${CMAKE_CURRENT_BINARY_DIR}/visualization.py COPYONLY)
configure_file(08-visualization.ipynb ${CMAKE_CURRENT_BINARY_DIR}/08-visualization.ipynb COPYONLY)

add_custom_target(tutorials_08)

html_export(TUTORIAL tutorials_08 FILE "08-visualization.ipynb")
python_export(TUTORIAL tutorials_08 FILE "08-visualization.ipynb")
