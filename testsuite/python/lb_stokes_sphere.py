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
import espressomd
import espressomd.lbboundaries
import espressomd.shapes
import unittest as ut
import unittest_decorators as utx
import numpy as np

# Define the LB parameters
TIME_STEP = 0.5
AGRID = 0.6
KVISC = 6
DENS = 2.3
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': KVISC,
             'tau': TIME_STEP}
# System setup
radius = 7 * AGRID
box_width = 46 * AGRID
real_width = box_width + 2 * AGRID
box_length = 36 * AGRID
c_s = np.sqrt(1. / 3. * AGRID**2 / TIME_STEP**2)
v = [0, 0, 0.1 * c_s]  # The boundary slip


class Stokes:
    """
    Measure the force on a single sphere immersed in a fluid with fixed
    velocity boundary conditions created by four walls at finite distance.
    The force is compared to th analytical result F=6 pi eta r v
    i.e. the stokes force on the particles.

    We create a box of size box_width x box_width x box_length and
    place an object in the center. We measure the drag force
    in z direction. We create walls in the xz and yz plane at the box
    boundaries, where the velocity is fixed to v.
    """
    system = espressomd.System(box_l=[real_width, real_width, box_length])
    system.box_l = [real_width, real_width, box_length]
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=1.0)

    def tearDown(self):
        self.system.actors.clear()
        self.system.lbboundaries.clear()
        self.system.thermostat.turn_off()

    def test_stokes(self):
        # Setup walls
        walls = [None] * 4
        walls[0] = espressomd.lbboundaries.LBBoundary(shape=espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(1 + box_width)), velocity=v)
        walls[1] = espressomd.lbboundaries.LBBoundary(shape=espressomd.shapes.Wall(
            normal=[1, 0, 0], dist=1), velocity=v)
        walls[2] = espressomd.lbboundaries.LBBoundary(shape=espressomd.shapes.Wall(
            normal=[0, -1, 0], dist=-(1 + box_width)), velocity=v)
        walls[3] = espressomd.lbboundaries.LBBoundary(shape=espressomd.shapes.Wall(
            normal=[0, 1, 0], dist=1), velocity=v)

        for wall in walls:
            self.system.lbboundaries.add(wall)

        # setup sphere without slip in the middle
        sphere = espressomd.lbboundaries.LBBoundary(shape=espressomd.shapes.Sphere(
            radius=radius, center=[real_width / 2] * 2 + [box_length / 2],
            direction=1))

        self.system.lbboundaries.add(sphere)

        def size(vector):
            tmp = 0
            for k in vector:
                tmp += k * k
            return np.sqrt(tmp)

        last_force = -1000.
        dynamic_viscosity = self.lbf.viscosity * self.lbf.density
        stokes_force = 6 * np.pi * dynamic_viscosity * radius * size(v)
        self.system.integrator.run(50)
        while True:
            self.system.integrator.run(3)
            force = np.linalg.norm(sphere.get_force())
            if np.abs(last_force - force) < 0.01 * stokes_force:
                break
            last_force = force

        force = np.copy(sphere.get_force())
        np.testing.assert_allclose(
            force,
            [0, 0, stokes_force],
            rtol=0.03,
            atol=stokes_force * 0.03)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(['LB_BOUNDARIES_GPU', 'EXTERNAL_FORCES'])
class LBGPUStokes(Stokes, ut.TestCase):
    lb_class = espressomd.lb.LBFluidGPU


@utx.skipIfMissingFeatures(['LB_BOUNDARIES', 'EXTERNAL_FORCES'])
class LBCPUStokes(Stokes, ut.TestCase):
    lb_class = espressomd.lb.LBFluid


if __name__ == "__main__":
    ut.main()
