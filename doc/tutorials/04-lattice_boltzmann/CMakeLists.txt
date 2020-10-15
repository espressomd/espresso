configure_tutorial_target(
  TARGET tutorial_04 DEPENDS 04-lattice_boltzmann_part1.ipynb
  04-lattice_boltzmann_part2.ipynb 04-lattice_boltzmann_part3.ipynb
  04-lattice_boltzmann_part4.ipynb figures/latticeboltzmann-grid.png
  figures/latticeboltzmann-momentumexchange.png)

nb_export(TARGET tutorial_04 SUFFIX "1" FILE "04-lattice_boltzmann_part1.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_04 SUFFIX "2" FILE "04-lattice_boltzmann_part2.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_04 SUFFIX "3" FILE "04-lattice_boltzmann_part3.ipynb"
          HTML_RUN VAR_SUBST "STEPS=200;BOX_L=16;N_MONOMERS=[6,8,10,12,14]")
nb_export(TARGET tutorial_04 SUFFIX "4" FILE "04-lattice_boltzmann_part4.ipynb"
          HTML_RUN)
