#
# Copyright (C) 2019-2022 The ESPResSo project
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
import sys
import numpy as np
import espressomd

system = espressomd.System(box_l=[100, 100, 100])
system.time_step = 0.01
system.cell_system.skin = 0.1

for i in range(100):
    system.part.add(pos=np.random.random() * system.box_l)

print("start of integration loop", file=sys.stderr)
while True:
    system.integrator.run(100000)
