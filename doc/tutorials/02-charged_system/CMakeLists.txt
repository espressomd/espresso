configure_file(scripts/nacl.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/nacl.py COPYONLY)
configure_file(scripts/nacl_units_confined.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/nacl_units_confined.py COPYONLY)
configure_file(scripts/nacl_units_confined_vis.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/nacl_units_confined_vis.py COPYONLY)
configure_file(scripts/nacl_units.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/nacl_units.py COPYONLY)
configure_file(scripts/nacl_units_vis.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/nacl_units_vis.py COPYONLY)
configure_file(figures/nacl-rdf.pdf ${CMAKE_CURRENT_BINARY_DIR}/figures/nacl-rdf.pdf COPYONLY)
configure_file(figures/nacl_units_confined.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/nacl_units_confined.jpg COPYONLY)
configure_file(figures/nacl_units.jpg ${CMAKE_CURRENT_BINARY_DIR}/figures/nacl_units.jpg COPYONLY)
configure_file(figures/salt.png ${CMAKE_CURRENT_BINARY_DIR}/figures/salt.png COPYONLY)
configure_file(02-charged_system-1.ipynb ${CMAKE_CURRENT_BINARY_DIR}/02-charged_system-1.ipynb COPYONLY)
configure_file(02-charged_system-2.ipynb ${CMAKE_CURRENT_BINARY_DIR}/02-charged_system-2.ipynb COPYONLY)

add_custom_target(tutorials_02)

nb_export(TUTORIAL tutorials_02_1 FILE "02-charged_system-1.ipynb" HTML_RUN
          VAR_SUBST "num_steps_equilibration=200;num_configs=100;integ_steps_per_config=250")
nb_export(TUTORIAL tutorials_02_2 FILE "02-charged_system-2.ipynb" HTML_RUN
          VAR_SUBST "num_steps_equilibration=200;num_configs=50;integ_steps_per_config=100")
