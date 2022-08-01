configure_tutorial_target(TARGET tutorial_widom DEPENDS widom_insertion.ipynb)

nb_export(TARGET tutorial_widom SUFFIX "" FILE "widom_insertion.ipynb" HTML_RUN)
