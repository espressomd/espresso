#
# Copyright (C) 2016-2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

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
