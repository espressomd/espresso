configure_file(05-raspberry_electrophoresis.ipynb ${CMAKE_CURRENT_BINARY_DIR}/05-raspberry_electrophoresis.ipynb COPYONLY)
configure_file(figures/raspberry_snapshot.png ${CMAKE_CURRENT_BINARY_DIR}/figures/raspberry_snapshot.png COPYONLY)

add_custom_target(tutorials_05)

html_export(TUTORIAL tutorials_05 FILE 05-raspberry_electrophoresis.ipynb)
python_export(TUTORIAL tutorials_05 FILE 05-raspberry_electrophoresis.ipynb)
