configure_tutorial_target(TARGET tutorial_08 DEPENDS 08-visualization.ipynb)

nb_export(TARGET tutorial_08 SUFFIX "" FILE "08-visualization.ipynb" HTML_RUN)
