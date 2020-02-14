configure_file(05-raspberry_electrophoresis.ipynb ${CMAKE_CURRENT_BINARY_DIR}/05-raspberry_electrophoresis.ipynb COPYONLY)
configure_file(figures/raspberry_snapshot.png ${CMAKE_CURRENT_BINARY_DIR}/figures/raspberry_snapshot.png COPYONLY)
configure_file(posVsTime_sample.dat ${CMAKE_CURRENT_BINARY_DIR}/posVsTime_sample.dat COPYONLY)

add_custom_target(tutorial_05)
add_dependencies(tutorials tutorial_05)

nb_export(TARGET tutorial_05 SUFFIX "" FILE 05-raspberry_electrophoresis.ipynb
          HTML_RUN VAR_SUBST "num_iterations=20;num_steps_per_iteration=20;'trajectory_file=\"posVsTime_sample.dat\"'")
