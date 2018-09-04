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
# Measuring the force on a single sphere immersed in a fluid with
# fixed velocity boundary conditions created by two
# walls at finite distance.
# The force is compared to th analytical result F=6 pi eta r v
# i.e. the stokes force on the particles.


# We create a box of size box_width x box_width x box_length and
# place an object in the center. We measure the drag force
# in z direction. We create walls in the xz and yz plane at the box
# boundaries, where the velocity is fixed to $v.
#
import espressomd
from espressomd import lb, lbboundaries, shapes, has_features
import unittest as ut
import numpy as np
import sys


@ut.skipIf(not has_features(["LB_GPU", "LB_BOUNDARIES_GPU"]),
           "Features not available, skipping test!")
class Stokes(ut.TestCase):

    def test_stokes(self):
        # System setup
        agrid = 1
        radius = 5.5
        box_width = 54
        real_width = box_width + 2 * agrid
        box_length = 54
        system = espressomd.System(box_l=[real_width, real_width, box_length])
        system.box_l = [real_width, real_width, box_length]
        system.time_step = 0.4
        system.cell_system.skin = 0.4

        # The temperature is zero.
        system.thermostat.set_lb(kT=0)

        # LB Parameters
        v = [0, 0, 0.01]  # The boundary slip
        kinematic_visc = 5.0

        # Invoke LB fluid
        lbf = lb.LBFluidGPU(visc=kinematic_visc, dens=1,
                            agrid=agrid, tau=system.time_step, fric=1)
        system.actors.add(lbf)

        # Setup walls
        walls = [None] * 4
        walls[0] = lbboundaries.LBBoundary(shape=shapes.Wall(
            normal=[-1, 0, 0], dist=-(1 + box_width)), velocity=v)
        walls[1] = lbboundaries.LBBoundary(
            shape=shapes.Wall(
                normal=[
                    1,
                    0,
                    0],
                dist=1),
            velocity=v)
        walls[2] = lbboundaries.LBBoundary(shape=shapes.Wall(
            normal=[0, -1, 0], dist=-(1 + box_width)), velocity=v)
        walls[3] = lbboundaries.LBBoundary(
            shape=shapes.Wall(
                normal=[
                    0,
                    1,
                    0],
                dist=1),
            velocity=v)

        for wall in walls:
            system.lbboundaries.add(wall)

        # setup sphere without slip in the middle
        sphere = lbboundaries.LBBoundary(shape=shapes.Sphere(
            radius=radius, center=[real_width / 2] * 2 + [box_length / 2], direction=1))

        system.lbboundaries.add(sphere)

        def size(vector):
            tmp = 0
            for k in vector:
                tmp += k * k
            return np.sqrt(tmp)

        system.integrator.run(800)

        stokes_force = 6 * np.pi * kinematic_visc * radius * size(v)
        print("Stokes' Law says: f=%f" % stokes_force)
        
        # get force that is exerted on the sphere
        for i in range(5):
            system.integrator.run(200)
            force = sphere.get_force()
            print("Measured force: f=%f" % size(force))
            self.assertLess(abs(1.0 - size(force) / stokes_force), 0.06)


if __name__ == "__main__":
    ut.main()
