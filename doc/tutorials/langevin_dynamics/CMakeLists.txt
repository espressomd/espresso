configure_tutorial_target(TARGET tutorial_ld DEPENDS langevin_dynamics.ipynb)

nb_export(TARGET tutorial_ld SUFFIX "" FILE "langevin_dynamics.ipynb" HTML_RUN)
