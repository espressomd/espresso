configure_tutorial_target(TARGET tutorial_01 DEPENDS 01-lennard_jones.ipynb
                          figures/lennard-jones-potential.png)

nb_export(TARGET tutorial_01 SUFFIX "" FILE "01-lennard_jones.ipynb" HTML_RUN
          VAR_SUBST)
