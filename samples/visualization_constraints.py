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
""" Visualization of shape-based constraints with test particles.
"""

from __future__ import print_function
from espressomd import *
from espressomd.shapes import *
from espressomd.visualization_opengl import *
from threading import Thread
import numpy as np

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

box_l = 50
system = espressomd.System(box_l=[50.0] * 3)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)

system.time_step = 0.0001
system.cell_system.skin = 0.3

visualizer = openGLLive(
    system, background_color=[1, 1, 1], drag_enabled=True, rasterize_resolution=50.0,
                        rasterize_pointsize=5, camera_position=[150, 25, 25], camera_right=[0, 0, -1])

# Wall
# system.constraints.add(shape=Wall(dist=20, normal=[0.1, 0.0, 1]),
# particle_type=0, penetrable=True)

# Sphere
# system.constraints.add(shape=Sphere(center=[25, 25, 25], radius=15,
# direction=1), particle_type=0, penetrable=True)

# Ellipsoid
# system.constraints.add(shape=Ellipsoid(center=[25, 25, 25], a=25, b=15,
# direction=1), particle_type=0, penetrable=True)

# Cylinder
# system.constraints.add(shape=Cylinder(center=[25, 25, 25], axis=[1, 0,
# 0], direction=1, radius=10, length=30), particle_type=0,
# penetrable=True)

# SpheroCylinder
# system.constraints.add(shape=SpheroCylinder(center=[25, 25, 25],
# axis=[1, 0, 0], direction=1, radius=10, length=30), particle_type=0,
# penetrable=True)

# Maze
# system.constraints.add(shape=Maze(cylrad=3, dim=2, nsphere=2, sphrad=8),
# particle_type=0, penetrable=True)

# Stomatocyte
# system.constraints.add(shape=Stomatocyte(inner_radius=3, outer_radius=7,
# axis=[1.0, 0.0, 0.0], center=[25, 25, 25], layer_width=3, direction=1),
# particle_type=0, penetrable=True)

# SimplePore
system.constraints.add(shape=SimplePore(center=[25, 25, 25], axis=[
                       1, 0, 0], length=15, radius=12.5, smoothing_radius=2), particle_type=0, penetrable=True)

# Slitpore
# system.constraints.add(shape=Slitpore(channel_width=15,
# lower_smoothing_radius=3, upper_smoothing_radius=3, pore_length=20,
# pore_mouth=30, pore_width=5), particle_type=0, penetrable=True)

# HollowCone
# system.constraints.add(shape=HollowCone(inner_radius=5, outer_radius=20,
# opening_angle=np.pi/4.0, axis=[1.0, 0.0, 0.0], center=[25, 25, 25],
# width=2, direction=1), particle_type=0, penetrable=True)


system.thermostat.set_langevin(kT=10.0, gamma=10)

for i in range(100):
    rpos = np.random.random(3) * box_l
    system.part.add(pos=rpos, type=1)

system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=1.0, sigma=5.0,
    cutoff=15.0, shift="auto")

system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=200.0, sigma=5.0,
    cutoff=20.0, shift="auto")

system.force_cap = 1000.0

visualizer.run(1)
