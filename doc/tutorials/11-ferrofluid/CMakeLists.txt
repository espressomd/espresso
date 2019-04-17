configure_file(11-ferrofluid_monolayer_part1.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_monolayer_part1.ipynb COPYONLY)
configure_file(11-ferrofluid_monolayer_part2.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_monolayer_part2.ipynb COPYONLY)
configure_file(11-ferrofluid_monolayer_part3.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_monolayer_part3.ipynb COPYONLY)
configure_file(figures/electrostatically_stabilised_2.svg ${CMAKE_CURRENT_BINARY_DIR}/figures/electrostatically_stabilised_2.svg COPYONLY)
configure_file(figures/Ferrofluid_Magnet_under_glass_edit.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/Ferrofluid_Magnet_under_glass_edit.jpg COPYONLY)
configure_file(figures/ferrofluid_system_schematic_Odenbach_Collodial_Magnetic_Fluids_page_261.png ${CMAKE_CURRENT_BINARY_DIR}/figures/ferrofluid_system_schematic_Odenbach_Collodial_Magnetic_Fluids_page_261.png COPYONLY)
configure_file(figures/headtotailconf.png ${CMAKE_CURRENT_BINARY_DIR}/figures/headtotailconf.png COPYONLY)
configure_file(figures/schematic_repr_sterically_stabilized.png ${CMAKE_CURRENT_BINARY_DIR}/figures/schematic_repr_sterically_stabilized.png COPYONLY)

add_custom_target(tutorials_11)

html_export(TUTORIAL tutorials_11_1 FILE "11-ferrofluid_monolayer_part1.ipynb")
html_export(TUTORIAL tutorials_11_2 FILE "11-ferrofluid_monolayer_part2.ipynb")
html_export(TUTORIAL tutorials_11_3 FILE "11-ferrofluid_monolayer_part3.ipynb")

python_export(TUTORIAL tutorials_11_1 FILE "11-ferrofluid_monolayer_part1.ipynb")
python_export(TUTORIAL tutorials_11_2 FILE "11-ferrofluid_monolayer_part2.ipynb")
python_export(TUTORIAL tutorials_11_3 FILE "11-ferrofluid_monolayer_part3.ipynb")
