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

import espressomd
import espressomd.shapes
import unittest as ut
import unittest_decorators as utx
import numpy as np

# Define the LB Parameters
TIME_STEP = 0.01
AGRID = 0.5
KVISC = 4
DENS = 2
G = 0.08
BOX_SIZE = 18 * AGRID

LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'viscosity': KVISC,
             'tau': TIME_STEP,
             'ext_force_density': [0, DENS * G, 0]}
# System setup
RADIUS = 6 * AGRID


class LBBuoyancy:
    """
    Tests buoyancy force on a sphere in a closed box of lb fluid and
    the overall force balance

    """
    system = espressomd.System(box_l=[BOX_SIZE] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()

    def test(self):
        # Setup walls
        for i in range(3):
            n = np.zeros(3)
            n[i] = 1
            self.lbf.add_boundary_from_shape(espressomd.shapes.Wall(
                normal=-n, dist=-(self.system.box_l[i] - AGRID)))
            self.lbf.add_boundary_from_shape(
                espressomd.shapes.Wall(normal=n, dist=AGRID))

        # setup sphere without slip in the middle
        sphere_shape = espressomd.shapes.Sphere(
            radius=RADIUS, center=self.system.box_l / 2, direction=1)
        self.lbf.add_boundary_from_shape(sphere_shape)

        sphere_volume = 4. / 3. * np.pi * RADIUS**3

        # TODO WALBERLA: (#4381)
        self.skipTest("boundary forces not implemented at the moment")

        # Equilibration
        last_force = np.inf * np.ones(3)
        self.system.integrator.run(100)
        while True:
            self.system.integrator.run(10)
            force = np.linalg.norm(self.lbf.boundary['sphere'].get_force())
            if np.linalg.norm(force - last_force) < 0.01:
                break
            last_force = force

        # Check force balance
        boundary_force = np.zeros(3)
        for b in self.lbf.boundary:
            boundary_force += b.get_force()

        fluid_nodes = np.sum(np.logical_not(
            self.lbf[:, :, :].is_boundary).astype(int))
        fluid_volume = fluid_nodes * AGRID**3
        applied_force = fluid_volume * np.array(LB_PARAMS['ext_force_density'])

        np.testing.assert_allclose(
            boundary_force,
            applied_force,
            atol=0.08 * np.linalg.norm(applied_force))

        # Check buoyancy force on the sphere
        expected_force = np.array(
            [0, -sphere_volume * DENS * G, 0])
        np.testing.assert_allclose(
            np.copy(self.lbf.boundary['sphere'].get_force()), expected_force,
            atol=np.linalg.norm(expected_force) * 0.02)


@utx.skipIfMissingFeatures(["EXTERNAL_FORCES", "WALBERLA"])
class LBBuoyancyWalberla(LBBuoyancy, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(["EXTERNAL_FORCES", "WALBERLA"])
class LBBuoyancyWalberlaSinglePrecision(LBBuoyancy, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == "__main__":
    ut.main()
