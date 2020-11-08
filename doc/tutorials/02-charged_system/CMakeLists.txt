configure_tutorial_target(
  TARGET
  tutorial_02
  DEPENDS
  02-charged_system-1.ipynb)

nb_export(TARGET tutorial_02 SUFFIX "1" FILE "02-charged_system-1.ipynb"
          HTML_RUN)
