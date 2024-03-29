#
# Copyright (C) 2010-2022 The ESPResSo project
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
Visualize the simulation of a linear polymer.
"""

import espressomd
import espressomd.interactions
import espressomd.visualization
import numpy as np

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

box_l = 50
n_part = 200

system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=0, sigma=1, cutoff=2, shift="auto")
system.bonded_inter[0] = espressomd.interactions.HarmonicBond(k=0.5, r_0=1.0)

previous_part = None
for i in range(n_part):
    part = system.part.add(pos=np.random.random(3) * system.box_l)
    if previous_part:
        part.add_bond((system.bonded_inter[0], previous_part))
    previous_part = part

visualizer = espressomd.visualization.openGLLive(
    system, bond_type_radius=[0.3])

system.integrator.set_steepest_descent(f_max=10, gamma=50.0,
                                       max_displacement=0.2)
system.integrator.run(1000)
system.integrator.set_vv()

system.thermostat.set_langevin(kT=0.1, gamma=20.0, seed=42)

visualizer.run(1)
