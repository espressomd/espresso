#
# Copyright (C) 2013-2019 The ESPResSo project
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
"""
Set up a lattice-Boltzmann fluid and apply an external force density on it.
"""
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(epilog=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument('--cpu', action='store_true')
group.add_argument('--gpu', action='store_true')
args = parser.parse_args()


print("""
=======================================================
=         lattice-Boltzmann fluid example             =
=======================================================
""")

required_features = ["EXTERNAL_FORCES"]
if args.gpu:
    print("Using GPU implementation")
    required_features.append("CUDA")
else:
    print("Using CPU implementation")
    if not args.cpu:
        print("(select the implementation with --cpu or --gpu)")

import espressomd
espressomd.assert_features(required_features)
import espressomd.lb


box_l = 50
system = espressomd.System(box_l=[box_l] * 3)

system.time_step = 0.01
system.cell_system.skin = 0.1

system.part.add(pos=[box_l / 2.0] * 3, fix=[True, True, True])


lb_params = {'agrid': 1, 'dens': 1, 'visc': 1, 'tau': 0.01,
             'ext_force_density': [0, 0, -1.0 / (box_l**3)]}

if args.gpu:
    lbf = espressomd.lb.LBFluidGPU(**lb_params)
else:
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
