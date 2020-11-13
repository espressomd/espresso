configure_tutorial_target(TARGET tutorial_viz DEPENDS visualization.ipynb)

nb_export(TARGET tutorial_viz SUFFIX "" FILE "visualization.ipynb" HTML_RUN)
