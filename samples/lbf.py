"""
This samples sets up a Lattice-Boltzmann fluid and applies an external force density on it.
"""

#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np

import espressomd

required_features = ["LB"]
espressomd.assert_features(required_features)

import espressomd.lb


print("""
=======================================================
=         Lattice Boltzmann fluid example             =
=======================================================

Program Information:""")
print(espressomd.features())


box_l = 50
system = espressomd.System(box_l=[box_l] * 3)
system.set_random_state_PRNG()

system.time_step = 0.01
system.cell_system.skin = 0.1

system.part.add(pos=[box_l / 2.0] * 3, fix=[1, 1, 1])


lb_params = {'agrid': 1, 'dens': 1, 'visc': 1, 'tau': 0.01,
             'ext_force_density': [0, 0, -1.0 / (box_l**3)]}
#lbf = espressomd.lb.LBFluidGPU(**lb_params)
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.0)
print(lbf.get_params())

f_list = np.zeros((10, 3))
for i in range(10):
    f_list[i] = system.part[0].f
    system.integrator.run(steps=10)
    print(i)

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.plot(f_list[:, 0], label=r"$F_x$")
ax.plot(f_list[:, 1], label=r"$F_y$")
ax.plot(f_list[:, 2], label=r"$F_z$")
ax.legend()
ax.set_xlabel("t")
ax.set_ylabel("F")

plt.show()
