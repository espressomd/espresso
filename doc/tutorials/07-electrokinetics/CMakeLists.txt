configure_file(scripts/eof_analytical.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/eof_analytical.py COPYONLY)
configure_file(scripts/plot.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/plot.py COPYONLY)
configure_file(figures/profiles.png ${CMAKE_CURRENT_BINARY_DIR}/figures/profiles.png  COPYONLY)
configure_file(figures/schlitzpore_3d.png ${CMAKE_CURRENT_BINARY_DIR}/figures/schlitzpore_3d.png COPYONLY)
configure_file(07-electrokinetics.ipynb ${CMAKE_CURRENT_BINARY_DIR}/07-electrokinetics.ipynb COPYONLY)

add_custom_target(tutorials_07)

nb_export(TUTORIAL tutorials_07 FILE "07-electrokinetics.ipynb" HTML_RUN VAR_SUBST "integration_length=600;dt=0.5" ADD_SCRIPTS "scripts/eof_analytical.py;scripts/plot.py")
