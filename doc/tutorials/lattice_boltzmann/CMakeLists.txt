configure_tutorial_target(
  TARGET tutorial_lb DEPENDS lattice_boltzmann_theory.ipynb
  lattice_boltzmann_poiseuille_flow.ipynb lattice_boltzmann_sedimentation.ipynb
  figures/latticeboltzmann-grid.png
  figures/latticeboltzmann-momentumexchange.png)

nb_export(TARGET tutorial_lb SUFFIX "1" FILE "lattice_boltzmann_theory.ipynb"
          HTML_RUN)
nb_export(TARGET tutorial_lb SUFFIX "2" FILE
          "lattice_boltzmann_poiseuille_flow.ipynb" HTML_RUN)
nb_export(TARGET tutorial_lb SUFFIX "3" FILE
          "lattice_boltzmann_sedimentation.ipynb" HTML_RUN)
