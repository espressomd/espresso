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
Visualize shape-based constraints interacting with a Lennard-Jones gas.
"""

import numpy as np
import argparse

import espressomd
import espressomd.shapes
import espressomd.visualization_opengl

parser = argparse.ArgumentParser(epilog=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument("--wall", action="store_const", dest="shape", const="Wall",
                   default="Wall")
for shape in ("Sphere", "Ellipsoid", "Cylinder", "SpheroCylinder",
              "Stomatocyte", "SimplePore", "Slitpore", "HollowConicalFrustum"):
    group.add_argument("--" + shape.lower(), action="store_const",
                       dest="shape", const=shape)
args = parser.parse_args()


required_features = ["LENNARD_JONES"]
if args.shape == "Stomatocyte":
    required_features.append("EXPERIMENTAL_FEATURES")
espressomd.assert_features(required_features)

box_l = 50.0
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.0001
system.cell_system.skin = 0.3

visualizer = espressomd.visualization_opengl.openGLLive(
    system,
    background_color=[1, 1, 1],
    drag_enabled=True,
    quality_constraints=200,
    rasterize_resolution=50.0,
    rasterize_pointsize=5,
    camera_position=[150, 25, 25],
    camera_right=[0, 0, -1])

if args.shape == "Wall":
    system.constraints.add(shape=espressomd.shapes.Wall(
        dist=20, normal=[0.1, 0.0, 1]),
        particle_type=0, penetrable=True)

elif args.shape == "Sphere":
    system.constraints.add(shape=espressomd.shapes.Sphere(
        center=[25, 25, 25], radius=15, direction=1),
        particle_type=0, penetrable=True)

elif args.shape == "Ellipsoid":
    system.constraints.add(shape=espressomd.shapes.Ellipsoid(
        center=[25, 25, 25], a=25, b=15, direction=1),
        particle_type=0, penetrable=True)

elif args.shape == "Cylinder":
    system.constraints.add(shape=espressomd.shapes.Cylinder(
        center=[25] * 3, axis=[1, 0, 0], direction=1, radius=10, length=30),
        particle_type=0,
        penetrable=True)

elif args.shape == "SpheroCylinder":
    system.constraints.add(
        shape=espressomd.shapes.SpheroCylinder(center=[25] * 3, axis=[1, 0, 0],
                                               direction=1, radius=10, length=30),
        particle_type=0,
        penetrable=True)

elif args.shape == "Stomatocyte":
    system.constraints.add(shape=espressomd.shapes.Stomatocyte(
        inner_radius=3, outer_radius=7, axis=[1.0, 0.0, 0.0], center=[25] * 3,
        layer_width=3, direction=1), particle_type=0, penetrable=True)

elif args.shape == "SimplePore":
    system.constraints.add(shape=espressomd.shapes.SimplePore(
        center=[25, 25, 25], axis=[1, 0, 0], length=15, radius=12.5,
        smoothing_radius=2), particle_type=0, penetrable=True)

elif args.shape == "Slitpore":
    system.constraints.add(shape=espressomd.shapes.Slitpore(
        channel_width=15, lower_smoothing_radius=3, upper_smoothing_radius=3,
        pore_length=20, pore_mouth=30, pore_width=5), particle_type=0,
        penetrable=True)

elif args.shape == "HollowConicalFrustum":
    system.constraints.add(shape=espressomd.shapes.HollowConicalFrustum(
        r1=12, r2=8, length=15.0, thickness=3,
        axis=[0.0, 1.0, 1.0], center=[25, 25, 25], direction=1),
        particle_type=0, penetrable=True)

else:
    raise ValueError("Unknown shape '{}'".format(args.shape))

for i in range(100):
    rpos = np.random.random(3) * box_l
    system.part.add(pos=rpos, type=1)

system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=1.0, sigma=5.0, cutoff=15.0, shift="auto")

system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=20.0, sigma=5.0, cutoff=20.0, shift="auto")


system.integrator.set_steepest_descent(f_max=10, gamma=1e-3,
                                       max_displacement=0.05)
system.integrator.run(500)
system.integrator.set_vv()
system.thermostat.set_langevin(kT=10.0, gamma=10, seed=42)

system.force_cap = 1000.0

visualizer.run(1)
