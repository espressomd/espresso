configure_file(11-ferrofluid_part1.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_part1.ipynb COPYONLY)
configure_file(11-ferrofluid_part2.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_part2.ipynb COPYONLY)
configure_file(11-ferrofluid_part3.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_part3.ipynb COPYONLY)
configure_file(figures/Electro-Steric_Stabilization.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/Electro-Steric_Stabilization.jpg COPYONLY)
configure_file(figures/Ferrofluid_Magnet_under_glass_edit.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/Ferrofluid_Magnet_under_glass_edit.jpg COPYONLY)
configure_file(figures/headtotailconf.png ${CMAKE_CURRENT_BINARY_DIR}/figures/headtotailconf.png COPYONLY)

add_custom_target(tutorials_11)

html_export(TUTORIAL tutorials_11_1 FILE "11-ferrofluid_part1.ipynb")
html_export(TUTORIAL tutorials_11_2 FILE "11-ferrofluid_part2.ipynb")
html_export(TUTORIAL tutorials_11_3 FILE "11-ferrofluid_part3.ipynb")

python_export(TUTORIAL tutorials_11_1 FILE "11-ferrofluid_part1.ipynb")
python_export(TUTORIAL tutorials_11_2 FILE "11-ferrofluid_part2.ipynb")
python_export(TUTORIAL tutorials_11_3 FILE "11-ferrofluid_part3.ipynb")
