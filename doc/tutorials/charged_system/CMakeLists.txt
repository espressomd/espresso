configure_tutorial_target(TARGET tutorial_crg DEPENDS charged_system.ipynb)

nb_export(TARGET tutorial_crg FILE "charged_system.ipynb" HTML_RUN)
