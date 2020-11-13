configure_tutorial_target(
  TARGET tutorial_lb DEPENDS lattice_boltzmann_part1.ipynb
  lattice_boltzmann_part2.ipynb lattice_boltzmann_part3.ipynb
  lattice_boltzmann_part4.ipynb figures/latticeboltzmann-grid.png
  figures/latticeboltzmann-momentumexchange.png)

nb_export(TARGET tutorial_lb SUFFIX "1" FILE "lattice_boltzmann_part1.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_lb SUFFIX "2" FILE "lattice_boltzmann_part2.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_lb SUFFIX "3" FILE "lattice_boltzmann_part3.ipynb"
          HTML_RUN VAR_SUBST "STEPS=200;BOX_L=16;N_MONOMERS=[6,8,10,12,14]")
nb_export(TARGET tutorial_lb SUFFIX "4" FILE "lattice_boltzmann_part4.ipynb"
          HTML_RUN)
