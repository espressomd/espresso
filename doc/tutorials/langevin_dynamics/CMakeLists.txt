configure_tutorial_target(TARGET langevin_dynamics DEPENDS
                          langevin_dynamics.ipynb polymers_part2.ipynb)

nb_export(TARGET tutorial_polymers SUFFIX "" FILE "langevin_dynamics.ipynb"
          HTML_RUN)
