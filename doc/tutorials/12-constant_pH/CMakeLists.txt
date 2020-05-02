configure_tutorial_target(TARGET tutorial_12 DEPENDS 12-constant_pH.ipynb)

nb_export(TARGET tutorial_12 SUFFIX "" FILE "12-constant_pH.ipynb" HTML_RUN)
