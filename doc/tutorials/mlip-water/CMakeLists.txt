#
# Copyright (C) 2026 The ESPResSo project
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

espresso_add_tutorial(
  TARGET tutorial_water SUFFIX 1 FILE "01_TIP4P_water.ipynb" DEPENDS
  soper2013tab5.csv TIP4P.png HTML_RUN VAR_SUBST
  "\"p3m_params={'cao':7,'mesh':[60,60,60]}\"")
espresso_add_tutorial(
  TARGET
  tutorial_water
  SUFFIX
  2
  FILE
  "02_MLIPs.ipynb"
  DEPENDS
  config/apax.yaml
  data.dvc
  dvc.yaml
  dvc.lock
  params.yaml
  pyproject.toml
  workflow.py)
espresso_add_tutorial(TARGET tutorial_water SUFFIX 3 FILE
                      "03_MLIP_within_ESPResSo.ipynb")
