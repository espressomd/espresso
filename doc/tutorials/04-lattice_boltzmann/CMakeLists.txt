configure_file(04-lattice_boltzmann_part4.ipynb ${CMAKE_CURRENT_BINARY_DIR}/04-lattice_boltzmann_part4.ipynb COPYONLY)
configure_file(04-lattice_boltzmann_part3.ipynb ${CMAKE_CURRENT_BINARY_DIR}/04-lattice_boltzmann_part3.ipynb COPYONLY)
configure_file(04-lattice_boltzmann_part2.ipynb ${CMAKE_CURRENT_BINARY_DIR}/04-lattice_boltzmann_part2.ipynb COPYONLY)
configure_file(04-lattice_boltzmann_part1.ipynb ${CMAKE_CURRENT_BINARY_DIR}/04-lattice_boltzmann_part1.ipynb COPYONLY)
configure_file(scripts/04-lattice_boltzmann_part3_solution.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part3_solution.py COPYONLY)
configure_file(scripts/04-lattice_boltzmann_part4_solution.py ${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution.py COPYONLY)
configure_file(figures/latticeboltzmann-momentumexchange.png ${CMAKE_CURRENT_BINARY_DIR}/figures/latticeboltzmann-momentumexchange.png COPYONLY)
configure_file(figures/latticeboltzmann-grid.png ${CMAKE_CURRENT_BINARY_DIR}/figures/latticeboltzmann-grid.png COPYONLY)

add_custom_command(
      OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution_cut.py"
      DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution.py"
      COMMAND
      "cp"
      "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution.py"
      "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution_cut.py"
      COMMAND
      "sed"
      "-n"
      "-i"
      "'/# Extract fluid velocity along the x-axis/,$$p;1i\\import matplotlib.pyplot as plt'"
      "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution_cut.py"
      )

add_custom_target(tutorials_04)

nb_export(TUTORIAL tutorials_04_1 FILE "04-lattice_boltzmann_part1.ipynb" HTML_RUN)
nb_export(TUTORIAL tutorials_04_2 FILE "04-lattice_boltzmann_part2.ipynb" HTML_RUN)
nb_export(TUTORIAL tutorials_04_3 FILE "04-lattice_boltzmann_part3.ipynb" HTML_RUN)
nb_export(TUTORIAL tutorials_04_4 FILE "04-lattice_boltzmann_part4.ipynb" HTML_RUN ADD_SCRIPTS "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution_cut.py")
