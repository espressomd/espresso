configure_tutorial_target(
  TARGET tutorial_05 DEPENDS 05-raspberry_electrophoresis.ipynb
  figures/raspberry_snapshot.png posVsTime_sample.dat)

nb_export(
  TARGET
  tutorial_05
  SUFFIX
  ""
  FILE
  05-raspberry_electrophoresis.ipynb
  HTML_RUN
  VAR_SUBST
  "num_iterations=20;num_steps_per_iteration=20;'trajectory_file=\"posVsTime_sample.dat\"'"
)
