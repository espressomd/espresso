configure_tutorial_target(
  TARGET
  tutorial_06
  DEPENDS
  06-active_matter.ipynb
  figures/friction.svg
  figures/pusher-puller.svg
  figures/enhanced.svg
  figures/geometry.svg
  figures/rectify.svg
  figures/rectification.svg
  figures/flow_field.svg
  exercises/enhanced_diffusion.py
  exercises/flow_field.py
  exercises/rectification_geometry.py
  exercises/rectification_simulation.py
  solutions/enhanced_diffusion.py
  solutions/flow_field.py
  solutions/rectification_geometry.py
  solutions/rectification_simulation.py
  solutions/solutions.txt)

nb_export(TARGET tutorial_06 SUFFIX "" FILE "06-active_matter.ipynb" HTML_RUN)
