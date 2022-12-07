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

import unittest as ut
import unittest_decorators as utx
import numpy as np
import itertools

import espressomd
import espressomd.shapes
import espressomd.lb

AGRID = 1.5
VISC = 4.2
DENS = 1.3
FRIC = 1.4
TAU = 0.2
BOX_L = 18.0
TIME_STEP = TAU
LB_PARAMETERS = {
    'agrid': AGRID,
    'viscosity': VISC,
    'density': DENS,
    'tau': TAU
}
V_BOUNDARY = 0.2


def velocity_profile(x):
    return V_BOUNDARY / (BOX_L - 2. * AGRID) * (x - AGRID)


class LBInterpolation:

    """
    Couette flow profile along x in z-direction. Check that velocity at shear
    plane next to the resting boundary is zero.
    """
    system = espressomd.System(box_l=[BOX_L] * 3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TIME_STEP

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMETERS, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()

    def set_boundaries(self, velocity):
        """Place boundaries *not* exactly on a LB node."""
        wall_shape1 = espressomd.shapes.Wall(
            normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(BOX_L - AGRID))
        self.lbf.add_boundary_from_shape(wall_shape1)
        self.lbf.add_boundary_from_shape(wall_shape2, velocity)

    def test_interpolated_velocity(self):
        """
        Check that the interpolated LB fluid velocity is zero between boundary
        node and first fluid node.
        """
        self.set_boundaries([0.0, 0.0, V_BOUNDARY])
        self.system.integrator.run(250)

        # Check interpolated vel at upper boundary. The node position is at
        # box_l[0]-agrid/2.
        np.testing.assert_allclose(
            np.copy(self.lbf.get_interpolated_velocity(
                [self.system.box_l[0] - AGRID / 2, 0, 0])),
            np.array([0, 0, V_BOUNDARY]))

        # Check interpolated velocity involving boundary and neighboring node.
        # The boundary node index is lbf.shape[0]-1, so -2 refers to the
        # node in front of the boundary.
        node_next_to_boundary = self.lbf[self.lbf.shape[0] - 2, 0, 0]
        # The midpoint between the boundary and that node is box_l - agrid.
        np.testing.assert_allclose(
            np.copy(self.lbf.get_interpolated_velocity([BOX_L - AGRID, 0, 0])),
            ([0, 0, V_BOUNDARY] + np.copy(node_next_to_boundary.velocity)) / 2.,
            atol=1e-6)

        # Bulk
        for pos in itertools.product(
                np.arange(1.5 * AGRID, BOX_L - 1.5 * AGRID, 0.5 * AGRID),
                np.arange(0.5 * AGRID, BOX_L, AGRID),
                np.arange(0.5 * AGRID, BOX_L, AGRID)):
            np.testing.assert_allclose(
                self.lbf.get_interpolated_velocity(pos)[2],
                velocity_profile(pos[0]), atol=5e-5)

    def test_mach_limit_check(self):
        """
        Assert that the Mach number check fires an exception.

        """
        max_vel = 1.1 * self.lbf.mach_limit() * AGRID / TAU
        vbb = espressomd.lb.VelocityBounceBack([0, 0, max_vel])
        error_msg = 'Slip velocity exceeds Mach 0.35'

        with self.assertRaisesRegex(ValueError, error_msg):
            self.lbf[0, 0, 0].boundary = vbb
        self.assertIsNone(self.lbf[0, 0, 0].boundary)

        with self.assertRaisesRegex(ValueError, error_msg):
            shape = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
            self.lbf.add_boundary_from_shape(shape, vbb.velocity)
        self.assertIsNone(self.lbf[0, 0, 0].boundary)


@utx.skipIfMissingFeatures(['WALBERLA'])
class LBInterpolationWalberla(LBInterpolation, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(['WALBERLA'])
class LBInterpolationWalberlaSinglePrecision(LBInterpolation, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == "__main__":
    ut.main()
