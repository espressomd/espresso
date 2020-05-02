configure_tutorial_target(
  TARGET tutorial_11 DEPENDS 11-ferrofluid_part1.ipynb
  11-ferrofluid_part2.ipynb 11-ferrofluid_part3.ipynb
  figures/Electro-Steric_Stabilization.jpg
  figures/Ferrofluid_Magnet_under_glass_edit.jpg figures/headtotailconf.png)

nb_export(TARGET tutorial_11 SUFFIX "1" FILE "11-ferrofluid_part1.ipynb"
          HTML_RUN VAR_SUBST "equil_steps=100;equil_rounds=10")
nb_export(TARGET tutorial_11 SUFFIX "2" FILE "11-ferrofluid_part2.ipynb"
          HTML_RUN VAR_SUBST
          "equil_steps=100;equil_rounds=10;alphas=[0,1,2,3,4.5,8];loops=100")
nb_export(
  TARGET tutorial_11 SUFFIX "3" FILE "11-ferrofluid_part3.ipynb" HTML_RUN
  VAR_SUBST
  "equil_steps=100;equil_rounds=10;alphas=[0.5,1,2,4,8];loops=200;loops_m=100")
