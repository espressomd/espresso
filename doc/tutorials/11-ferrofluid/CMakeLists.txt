configure_file(11-ferrofluid_part1.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_part1.ipynb COPYONLY)
configure_file(11-ferrofluid_part2.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_part2.ipynb COPYONLY)
configure_file(11-ferrofluid_part3.ipynb ${CMAKE_CURRENT_BINARY_DIR}/11-ferrofluid_part3.ipynb COPYONLY)
configure_file(figures/Electro-Steric_Stabilization.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/Electro-Steric_Stabilization.jpg COPYONLY)
configure_file(figures/Ferrofluid_Magnet_under_glass_edit.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/Ferrofluid_Magnet_under_glass_edit.jpg COPYONLY)
configure_file(figures/headtotailconf.png ${CMAKE_CURRENT_BINARY_DIR}/figures/headtotailconf.png COPYONLY)

add_custom_target(tutorials_11)

nb_export(TUTORIAL tutorials_11_1 FILE "11-ferrofluid_part1.ipynb" HTML_RUN VAR_SUBST "equil_steps=100;equil_rounds=10")
nb_export(TUTORIAL tutorials_11_2 FILE "11-ferrofluid_part2.ipynb" HTML_RUN VAR_SUBST "equil_steps=100;equil_rounds=10;alphas=[0,1,2,3,4.5,8];loops=100")
nb_export(TUTORIAL tutorials_11_3 FILE "11-ferrofluid_part3.ipynb" HTML_RUN VAR_SUBST "equil_steps=100;equil_rounds=10;alphas=[0.5,1,2,4,8];loops=200;loops_m=100")
