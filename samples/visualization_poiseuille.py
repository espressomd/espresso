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
Visualization sample for Poiseuille flow with Lattice Boltzmann.
"""

from __future__ import print_function
from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from threading import Thread
import espressomd.visualization_opengl

required_features = ["LB", "LB_BOUNDARIES", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

# System setup
box_l = 16
system = System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.cell_system.skin = 0.2

visualizer = espressomd.visualization_opengl.openGLLive(
    system,
    LB_draw_boundaries=True,
    LB_draw_velocity_plane=True,
    LB_plane_dist=8,
    LB_plane_axis=1,
    LB_vel_scale=1e2,
    LB_plane_ngrid=15,
    camera_position=[8, 16, 50],
    velocity_arrows=True,
    velocity_arrows_type_scale=[20.],
    velocity_arrows_type_radii=[0.1],
    velocity_arrows_type_colors=[[0, 1, 0]])

lbf = lb.LBFluid(kT=0, agrid=1.0, dens=1.0,
                 visc=1.0, tau=0.1, ext_force_density=[0, 0.003, 0])
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.5)

# Setup boundaries
walls = [lbboundaries.LBBoundary() for k in range(2)]
walls[0].set_params(shape=shapes.Wall(normal=[1, 0, 0], dist=1.5))
walls[1].set_params(shape=shapes.Wall(normal=[-1, 0, 0], dist=-14.5))

for i in range(100):
    system.part.add(pos=np.random.random(3) * system.box_l)

for wall in walls:
    system.lbboundaries.add(wall)

visualizer.run(1)
