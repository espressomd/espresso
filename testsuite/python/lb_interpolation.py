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
import unittest as ut
import unittest_decorators as utx
import numpy as np
import itertools
import sys

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
    'visc': VISC,
    'dens': DENS,
    'tau': TAU
}
V_BOUNDARY = 0.6


def velocity_profile(x):
    return V_BOUNDARY / (BOX_L - 2. * AGRID) * (x - AGRID)


class LBInterpolation:

    """
    Couette flow profile along x in z-direction. Check that velocity at shear
    plane next to the resting boundary is zero.
    """
    lbf = None
    system = espressomd.System(box_l=[BOX_L] * 3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TIME_STEP

    def set_boundaries(self, velocity):
        """Place boundaries *not* exactly on a LB node."""
        wall_shape1 = espressomd.shapes.Wall(
            normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(BOX_L - AGRID))
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=wall_shape1))
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=wall_shape2, velocity=velocity))

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

        # Check interpolated velocity involving boundary and neighboring node
        # The boundary node index is lbf.shape[0]-1, so -2 refers to the
        # node in front of the boundary. 
        node_next_to_boundary = self.lbf[
            self.lbf.shape[0] - 2, 0, 0]
        # The midpoint between the boundary and
        # that node is box_l - agrid.
        np.testing.assert_allclose(
            np.copy(self.lbf.get_interpolated_velocity(
                [self.system.box_l[0] - AGRID, 0, 0])),
            0.5 * (np.array([0, 0, V_BOUNDARY]) + node_next_to_boundary.velocity))

        # Bulk
        for pos in itertools.product(
                np.arange(1.5 * AGRID, BOX_L - 1.5 * AGRID, 0.5 * AGRID),
                np.arange(0.5 * AGRID, BOX_L, AGRID),
                np.arange(0.5 * AGRID, BOX_L, AGRID)):
            np.testing.assert_almost_equal(
                self.lbf.get_interpolated_velocity(pos)[2], velocity_profile(pos[0]), decimal=4)
        # Shear plane for boundary 2
        # for pos in itertools.product((9 * AGRID,), np.arange(0.5 * AGRID, BOX_L, AGRID), np.arange(0.5 * AGRID, BOX_L, AGRID)):
        # np.testing.assert_almost_equal(self.lbf.get_interpolated_velocity(pos)[2],
        # 1.0, decimal=4)

    def test_mach_limit_check(self):
        """
        Assert that the mach number check fires an exception.

        """
        max_vel = 0.31 * AGRID / TAU
        print("Begin: Test error generation")
        sys.stdout.flush()
        sys.stderr.flush()
        with self.assertRaises(Exception):
            self.set_boundaries([0.0, 0.0, max_vel])
            self.system.integrator.run(1)
        sys.stdout.flush()
        sys.stderr.flush()
        print("End: Test error generation")


@utx.skipIfMissingFeatures(['LB_BOUNDARIES'])
class LBInterpolationCPU(ut.TestCase, LBInterpolation):

    def setUp(self):
        self.system.lbboundaries.clear()
        self.system.actors.clear()
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMETERS)
        self.system.actors.add(self.lbf)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(['LB_BOUNDARIES_GPU'])
class LBInterpolationGPU(ut.TestCase, LBInterpolation):

    def setUp(self):
        self.system.lbboundaries.clear()
        self.system.actors.clear()
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMETERS)
        self.system.actors.add(self.lbf)


if __name__ == "__main__":
    ut.main()
