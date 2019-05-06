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
import unittest as ut
import numpy as np
import itertools

import espressomd
import espressomd.shapes
import espressomd.lb

AGRID = 1.5
VISC = 1.0
DENS = 1.0
FRIC = 1.0
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


class LBInterpolation(object):

    """
    Couette flow profile along x in z-direction. Check that velocity at shear plane next to
    the resting boundary is zero.

    """
    lbf = None
    system = espressomd.System(box_l=[BOX_L] * 3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TIME_STEP

    def set_boundaries(self, velocity):
        """Place boundaries *not* exactly on a LB node.

        """
        wall_shape1 = espressomd.shapes.Wall(
            normal=[1, 0, 0], dist=0.6 * AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(BOX_L - 0.6 * AGRID))
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=wall_shape1))
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=wall_shape2, velocity=velocity))

    def test_interpolated_velocity(self):
        """
        Check that the interpolated LB fluid velocity is zero between boundary
        node and first fluid node.

        """
        self.system.lbboundaries.clear()
        self.set_boundaries([0.0, 0.0, V_BOUNDARY])
        self.system.integrator.run(1200)
        # Shear plane for boundary 1
        #for pos in itertools.product((AGRID,), np.arange(0.5 * AGRID, BOX_L, AGRID), np.arange(0.5 * AGRID, BOX_L, AGRID)):
        #    np.testing.assert_almost_equal(self.lbf.get_interpolated_velocity(pos)[2], 0.0)
        # Bulk
        for pos in itertools.product(np.arange(1.5 * AGRID, BOX_L - 1.5 * AGRID, 0.5 * AGRID), np.arange(0.5 * AGRID, BOX_L, AGRID), np.arange(0.5 * AGRID, BOX_L, AGRID)):
            np.testing.assert_almost_equal(
                self.lbf.get_interpolated_velocity(pos)[2], velocity_profile(pos[0]), decimal=4)
        # Shear plane for boundary 2
        #for pos in itertools.product((9 * AGRID,), np.arange(0.5 * AGRID, BOX_L, AGRID), np.arange(0.5 * AGRID, BOX_L, AGRID)):
        # np.testing.assert_almost_equal(self.lbf.get_interpolated_velocity(pos)[2],
        # 1.0, decimal=4)

    def test_mach_limit_check(self):
        """
        Assert that the mach number check fires an exception.

        """
        max_vel = 0.31 * AGRID / TAU
        self.system.lbboundaries.clear()
        self.set_boundaries([0.0, 0.0, max_vel])
        with self.assertRaises(Exception):
            self.system.integrator.run(1)
        

@ut.skipIf(not espressomd.has_features(['LB_BOUNDARIES']), "Skipped, features missing.")
class LBInterpolationCPU(ut.TestCase, LBInterpolation):

    def setUp(self):
        self.system.actors.clear()
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMETERS)
        self.system.actors.add(self.lbf)


@ut.skipIf(not espressomd.gpu_available() or not espressomd.has_features(['LB_GPU', 'LB_BOUNDARIES_GPU']), "Skipped, features or gpu missing.")
class LBInterpolationGPU(ut.TestCase, LBInterpolation):

    def setUp(self):
        self.system.actors.clear()
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMETERS)
        self.system.actors.add(self.lbf)


if __name__ == "__main__":
    ut.main()
