configure_tutorial_target(TARGET tutorial_langevin_dynamics DEPENDS
                          langevin_dynamics.ipynb)

nb_export(TARGET tutorial_langevin_dynamics SUFFIX "" FILE
          "langevin_dynamics.ipynb" HTML_RUN)
