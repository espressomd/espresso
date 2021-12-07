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
"""
Stokesian Dynamics simulation of particle sedimentation.
Reproduce the trajectory in Figure 5b from :cite:`durlofsky87a`.
"""
import espressomd
import espressomd.constraints
import espressomd.observables
import espressomd.accumulators
import numpy as np
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(epilog=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument('--ft', action='store_true', help='Use FT approximation')
group.add_argument('--fts', action='store_true', help='Use FTS approximation')
args = parser.parse_args()

if args.ft:
    print("Using FT approximation method")
    sd_method = "ft"
else:
    print("Using FTS approximation method")
    sd_method = "fts"

espressomd.assert_features("STOKESIAN_DYNAMICS")

system = espressomd.System(box_l=[10, 10, 10])
system.time_step = 1.5
system.cell_system.skin = 0.4
system.periodicity = [False, False, False]

system.integrator.set_stokesian_dynamics(
    viscosity=1.0, radii={0: 1.0}, approximation_method=sd_method)

system.part.add(pos=[-5, 0, 0], rotation=[1, 1, 1])
system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])
system.part.add(pos=[7, 0, 0], rotation=[1, 1, 1])

gravity = espressomd.constraints.Gravity(g=[0, -1, 0])
system.constraints.add(gravity)

obs = espressomd.observables.ParticlePositions(ids=system.part.all().id)
acc = espressomd.accumulators.TimeSeries(obs=obs, delta_N=1)
system.auto_update_accumulators.add(acc)
acc.update()
intsteps = int(10500 / system.time_step)
system.integrator.run(intsteps)

positions = acc.time_series()
ref_data = "../testsuite/python/data/dancing.txt"
data = np.loadtxt(ref_data)

for i in range(3):
    plt.plot(positions[:, i, 0], positions[:, i, 1], linestyle='solid')

plt.gca().set_prop_cycle(None)

for i in range(0, 6, 2):
    plt.plot(data[:, i], data[:, i + 1], linestyle='dashed')

plt.title(f"Trajectory of sedimenting spheres\nsolid line: simulation "
          f"({sd_method.upper()}), dashed line: paper (FTS)")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
