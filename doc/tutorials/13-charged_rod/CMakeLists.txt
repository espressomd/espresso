configure_tutorial_target(TARGET tutorial_13 DEPENDS 13-charged_rod.ipynb)

nb_export(TARGET tutorial_13 SUFFIX "" FILE "13-charged_rod.ipynb" HTML_RUN)
