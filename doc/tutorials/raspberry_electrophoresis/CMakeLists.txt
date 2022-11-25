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
  TARGET tutorial_eph DEPENDS raspberry_electrophoresis.ipynb
  figures/raspberry_snapshot.png posVsTime_sample.dat)

nb_export(
  TARGET
  tutorial_eph
  SUFFIX
  ""
  FILE
  raspberry_electrophoresis.ipynb
  HTML_RUN
  VAR_SUBST
  "num_iterations=20;num_steps_per_iteration=20;'trajectory_file=\"posVsTime_sample.dat\"'"
)
