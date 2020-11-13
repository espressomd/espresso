configure_tutorial_target(TARGET tutorial_cph DEPENDS constant_pH.ipynb)

nb_export(TARGET tutorial_cph SUFFIX "" FILE "constant_pH.ipynb" HTML_RUN)
