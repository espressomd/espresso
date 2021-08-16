configure_tutorial_target(TARGET tutorial_polymers DEPENDS polymers.ipynb)

nb_export(TARGET tutorial_polymers SUFFIX "" FILE "polymers.ipynb" HTML_RUN
          VAR_SUBST "STEPS=200;BOX_L=16;N_MONOMERS=[6,8,10,12,14]")
