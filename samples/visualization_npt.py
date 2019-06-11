# Copyright (C) 2010-2018 The ESPResSo project
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
"""
Visualization sample for particle dumbbells in the constant-temperature,
constant-pressure ensemble.
"""

from __future__ import print_function
import numpy as np
from threading import Thread

import espressomd
from espressomd import thermostat
from espressomd.interactions import HarmonicBond
import espressomd.visualization_opengl

required_features = ["NPT", "LENNARD_JONES"]
espressomd.assert_features(required_features)

box_l = 10
system = espressomd.System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)

visualizer = espressomd.visualization_opengl.openGLLive(
    system, background_color=[1, 1, 1], bond_type_radius=[0.2])

system.time_step = 0.0005
system.cell_system.skin = 0.1

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=2, sigma=1,
    cutoff=3, shift="auto")

system.bonded_inter[0] = HarmonicBond(k=5.0, r_0=1.0)

n_part = 200
for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

for i in range(0, n_part - 1, 2):
    system.part[i].add_bond((system.bonded_inter[0], system.part[i + 1].id))

print("E before minimization:", system.analysis.energy()["total"])
system.minimize_energy.init(f_max=0.0, gamma=30.0,
                            max_steps=10000, max_displacement=0.1)
system.minimize_energy.minimize()
print("E after minimization:", system.analysis.energy()["total"])

system.thermostat.set_npt(kT=2.0, gamma0=1.0, gammav=0.01)
system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=0.01)


def main():
    cnt = 0
    P = 0
    while True:
        system.integrator.run(1)
        P += system.analysis.pressure()['total']
        if cnt > 10000:
            print("Pressure:", P / cnt, "Box:", system.box_l)
            cnt = 0
            P = 0

        visualizer.update()
        cnt += 1

# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
