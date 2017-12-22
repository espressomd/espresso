#
# Copyright (C) 2013,2014 The ESPResSo project
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
import espressomd.lb


print("""
=======================================================
=         Lattice Boltzmann fluid example             =
=======================================================

Program Information:""")
print(espressomd.features())


system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.1
box_l = 50
system.box_l = [box_l, box_l, box_l]

system.part.add(id=0, pos=[box_l / 2.0, box_l /
                           2.0, box_l / 2.0], fix=[1, 1, 1])


lb_params = {'agrid': 1, 'fric': 1, 'dens': 1, 'visc': 1, 'tau': 0.01,
             'ext_force': [0, 0, -1.0 / (box_l**3)]}
#lbf = espressomd.lb.LBFluidGPU(**lb_params)
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
print(system.actors)
print(lbf.get_params())

f_list = []
for i in range(10):
    f_list.append(system.part[0].f)
    system.integrator.run(steps=10)
    print(i)

f_list = np.array(f_list)


fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.plot(f_list[:, 0], label="F_x")
ax.plot(f_list[:, 1], label="F_y")
ax.plot(f_list[:, 2], label="F_z")
ax.legend()
ax.set_xlabel("t")
ax.set_ylabel("F")

plt.show()
