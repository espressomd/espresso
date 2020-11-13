configure_tutorial_target(TARGET tutorial_crg DEPENDS charged_system-1.ipynb)

nb_export(TARGET tutorial_crg SUFFIX "1" FILE "charged_system-1.ipynb" HTML_RUN)
