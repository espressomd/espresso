configure_tutorial_target(
  TARGET
  tutorial_02
  DEPENDS
  02-charged_system-1.ipynb
  02-charged_system-2.ipynb
  figures/nacl-rdf.pdf
  figures/nacl_units_confined.jpg
  figures/nacl_units.jpg
  figures/salt.png
  scripts/nacl.py
  scripts/nacl_units_confined.py
  scripts/nacl_units_confined_vis.py
  scripts/nacl_units.py
  scripts/nacl_units_vis.py)

nb_export(
  TARGET tutorial_02 SUFFIX "1" FILE "02-charged_system-1.ipynb" HTML_RUN
  VAR_SUBST
  "num_steps_equilibration=200;num_configs=100;integ_steps_per_config=250")
nb_export(
  TARGET tutorial_02 SUFFIX "2" FILE "02-charged_system-2.ipynb" HTML_RUN
  VAR_SUBST
  "num_steps_equilibration=200;num_configs=50;integ_steps_per_config=100")
