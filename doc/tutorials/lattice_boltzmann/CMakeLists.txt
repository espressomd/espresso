configure_tutorial_target(
  TARGET tutorial_lb DEPENDS lattice_boltzmann_part1.ipynb
  lattice_boltzmann_part2.ipynb figures/latticeboltzmann-grid.png
  figures/latticeboltzmann-momentumexchange.png)

nb_export(TARGET tutorial_lb SUFFIX "1" FILE "lattice_boltzmann_part1.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_lb SUFFIX "2" FILE "lattice_boltzmann_part2.ipynb"
          HTML_RUN)