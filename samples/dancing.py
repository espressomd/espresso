#
# Copyright (C) 2019-2020 The ESPResSo project
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

import espressomd
from espressomd import constraints
import numpy as np
import matplotlib.pyplot as plt

espressomd.assert_features("STOKESIAN_DYNAMICS")

system = espressomd.System(box_l=[10, 10, 10])
system.time_step = 1.0
system.cell_system.skin = 0.4

system.integrator.set_sd(viscosity=1.0, radii={0: 1.0})

system.part.add(pos=[-5, 0, 0], rotation=[1, 1, 1])
system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])
system.part.add(pos=[7, 0, 0], rotation=[1, 1, 1])

gravity = constraints.Gravity(g=[0, -1, 0])
system.constraints.add(gravity)

intsteps = int(13000 / system.time_step)
pos = np.empty([intsteps, 3 * len(system.part)])
for i in range(intsteps):
    system.integrator.run(1)
    for n, p in enumerate(system.part):
        pos[i, 3 * n:3 * n + 3] = p.pos

for n, p in enumerate(system.part):
    plt.plot(pos[:, 3 * n], pos[:, 3 * n + 1])
plt.show()
