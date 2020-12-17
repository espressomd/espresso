configure_tutorial_target(
  TARGET tutorial_polymers DEPENDS polymers_part1.ipynb
  polymers_part2.ipynb)

nb_export(TARGET tutorial_polymers SUFFIX "1" FILE "polymers_part1.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_polymers SUFFIX "3" FILE "polymers_part2.ipynb"
          HTML_RUN VAR_SUBST "STEPS=200;BOX_L=16;N_MONOMERS=[6,8,10,12,14]")