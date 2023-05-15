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
Visualize the Poiseuille flow in a lattice-Boltzmann fluid with an
external force applied.
"""

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.visualization
import numpy as np

required_features = ["WALBERLA", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

# System setup
box_l = 16
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.2

visualizer = espressomd.visualization.openGLLive(
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

lbf = espressomd.lb.LBFluidWalberla(kT=0, agrid=1.0, density=1.0, kinematic_viscosity=1.0,
                                    tau=0.1, ext_force_density=[0, 0.003, 0])
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.5)

# Setup boundaries
wall_shapes = [None] * 2
wall_shapes[0] = espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5)
wall_shapes[1] = espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-14.5)

for i in range(100):
    system.part.add(pos=np.random.random(3) * system.box_l)

for wall_shape in wall_shapes:
    lbf.add_boundary_from_shape(wall_shape)

visualizer.run(1)
