#
# Copyright (C) 2022-2023 The ESPResSo project
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

foreach(n_species RANGE 1 5)
  foreach(precision double_precision single_precision)
    target_sources(
      espresso_walberla
      PRIVATE ReactionKernelBulk_${n_species}_${precision}.cpp
              ReactionKernelIndexed_${n_species}_${precision}.cpp)
  endforeach()
endforeach()
