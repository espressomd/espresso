configure_tutorial_target(
  TARGET tutorial_06
  DEPENDS
  06-active_matter.ipynb
  FIGURES/friction.svg
  FIGURES/pusher-puller.svg
  FIGURES/enhanced.svg
  FIGURES/geometry.svg
  FIGURES/rectify.svg
  FIGURES/rectification.svg
  FIGURES/flow_field.svg
  EXERCISES/enhanced_diffusion.py
  EXERCISES/flow_field.py
  EXERCISES/rectification_geometry.py
  EXERCISES/rectification_simulation.py
  SOLUTIONS/enhanced_diffusion.py
  SOLUTIONS/flow_field.py
  SOLUTIONS/rectification_geometry.py
  SOLUTIONS/rectification_simulation.py
  SOLUTIONS/solutions.txt
)

nb_export(TARGET tutorial_06 SUFFIX "" FILE "06-active_matter.ipynb" HTML_RUN)
