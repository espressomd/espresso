configure_tutorial_target(
  TARGET tutorial_am DEPENDS active_matter.ipynb figures/friction.svg
  figures/pusher-puller.svg figures/geometry.svg figures/flow_field.svg)

nb_export(TARGET tutorial_am FILE "active_matter.ipynb" HTML_RUN)
