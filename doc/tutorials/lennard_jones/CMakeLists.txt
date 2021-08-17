configure_tutorial_target(TARGET tutorial_lj DEPENDS lennard_jones.ipynb)

nb_export(TARGET tutorial_lj SUFFIX "" FILE "lennard_jones.ipynb" HTML_RUN
          VAR_SUBST)
