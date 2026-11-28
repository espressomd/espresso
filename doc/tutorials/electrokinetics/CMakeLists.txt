#
# Copyright (C) 2016-2026 The ESPResSo project
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
  TARGET tutorial_ek SUFFIX "1" FILE
  "electrokinetics_part1_advection_diffusion.ipynb" HTML_RUN)
espresso_add_tutorial(
  TARGET tutorial_ek SUFFIX "2" FILE
  "electrokinetics_part2_electroosmotic_flow.ipynb" HTML_RUN DEPENDS
  figures/schlitzpore_3d.png)
espresso_add_tutorial(TARGET tutorial_ek SUFFIX "3" FILE
                      "electrokinetics_part3_reactive_flow.ipynb" HTML_RUN)
