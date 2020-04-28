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
Visualize lattice-Boltzmann boundary nodes.
"""

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.lbboundaries
from espressomd import visualization_opengl

required_features = ["LB_BOUNDARIES"]
espressomd.assert_features(required_features)

system = espressomd.System(box_l=[10.0, 10.0, 5.0])
system.time_step = 0.01
system.cell_system.skin = 0.4

lb_fluid = espressomd.lb.LBFluid(
    agrid=1.0, dens=1.0, visc=1.0, tau=0.01, ext_force_density=[0, 0, 0.15])
system.actors.add(lb_fluid)

cylinder_shape = espressomd.shapes.Cylinder(
    center=[5.0, 5.0, 5.0],
    axis=[0, 0, 1],
    direction=-1,
    radius=4.0,
    length=20.0)
cylinder_boundary = espressomd.lbboundaries.LBBoundary(shape=cylinder_shape)
system.lbboundaries.add(cylinder_boundary)

visualizer = visualization_opengl.openGLLive(
    system,
    background_color=[1, 1, 1],
    camera_position=[5, 5, 25],
    LB_draw_boundaries=True,
    LB_draw_nodes=True,
    LB_draw_node_boundaries=True)

visualizer.run(1)
