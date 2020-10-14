configure_tutorial_target(
  TARGET tutorial_04 DEPENDS 04-lattice_boltzmann_part1.ipynb
  04-lattice_boltzmann_part2.ipynb 04-lattice_boltzmann_part3.ipynb
  04-lattice_boltzmann_part4.ipynb figures/latticeboltzmann-grid.png
  figures/latticeboltzmann-momentumexchange.png
  scripts/04-lattice_boltzmann_part3_solution.py)

nb_export(TARGET tutorial_04 SUFFIX "1" FILE "04-lattice_boltzmann_part1.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_04 SUFFIX "2" FILE "04-lattice_boltzmann_part2.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_04 SUFFIX "3" FILE "04-lattice_boltzmann_part3.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_04 SUFFIX "4" FILE "04-lattice_boltzmann_part4.ipynb"
          HTML_RUN)
