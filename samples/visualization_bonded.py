# Copyright (C) 2010-2019 The ESPResSo project
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
Visualize the simulation of a linear polymer.
"""

import espressomd
from espressomd.interactions import HarmonicBond
from espressomd import visualization
from espressomd.minimize_energy import steepest_descent
import numpy as np
import argparse

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

parser = argparse.ArgumentParser(epilog=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument("--mayavi", action="store_const", dest="visualizer",
                   const="mayavi", help="MayaVi visualizer", default="mayavi")
group.add_argument("--opengl", action="store_const", dest="visualizer",
                   const="opengl", help="OpenGL visualizer")
args = parser.parse_args()

box_l = 50
n_part = 200

system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=0.1, gamma=20.0, seed=42)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=0, sigma=1, cutoff=2, shift="auto")
system.bonded_inter[0] = HarmonicBond(k=0.5, r_0=1.0)

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

for i in range(n_part - 1):
    system.part[i].add_bond((system.bonded_inter[0], system.part[i + 1].id))

# Select visualizer
if args.visualizer == "mayavi":
    visualizer = visualization.mayaviLive(system)
else:
    visualizer = visualization.openGLLive(system, bond_type_radius=[0.3])

steepest_descent(system, f_max=10, gamma=50.0, max_steps=1000,
                 max_displacement=0.2)

visualizer.run(1)
