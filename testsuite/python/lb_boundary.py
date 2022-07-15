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
import espressomd
import espressomd.lb
import espressomd.shapes
import numpy as np


class LBBoundariesBase:
    system = espressomd.System(box_l=[10.0, 5.0, 5.0])
    system.cell_system.skin = 0.1

    wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
    wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-7.5)

    def setUp(self):
        self.lbf = self.lb_class(
            viscosity=1.0, density=1.0, agrid=0.5, tau=1.0, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()

    def check_boundary_flags(self, slip_velocity1, slip_velocity2):
        lbb1 = self.lbf[:5, :, :]
        lbb2 = self.lbf[15:, :, :]
        lbb3 = self.lbf[5:15, :, :]
        np.testing.assert_equal(np.copy(lbb1.is_boundary), True)
        np.testing.assert_equal(np.copy(lbb2.is_boundary), True)
        np.testing.assert_equal(np.copy(lbb3.is_boundary), False)
        np.testing.assert_allclose(np.copy(lbb1.velocity),
                                   np.tile(slip_velocity1, [5, 10, 10, 1]))
        np.testing.assert_allclose(np.copy(lbb2.velocity),
                                   np.tile(slip_velocity2, [5, 10, 10, 1]))
        self.lbf.clear_boundaries()
        np.testing.assert_equal(np.copy(self.lbf[:, :, :].is_boundary), False)

    def test_boundary_flags(self):
        slip_velocity1 = 1e-3 * np.array([1., 2., 3.])
        slip_velocity2 = 1e-3 * np.array([4., 5., 6.])
        self.lbf.add_boundary_from_shape(self.wall_shape1, slip_velocity1)
        self.lbf.add_boundary_from_shape(self.wall_shape2, slip_velocity2)
        self.check_boundary_flags(slip_velocity1, slip_velocity2)

    def test_union(self):
        union = espressomd.shapes.Union()
        union.add([self.wall_shape1, self.wall_shape2])

        slip_velocity = 1e-3 * np.array([1., 2., 3.])
        self.lbf.add_boundary_from_shape(union, slip_velocity)
        self.check_boundary_flags(slip_velocity, slip_velocity)


@utx.skipIfMissingFeatures(["LB_WALBERLA"])
class LBBoundariesWalberla(LBBoundariesBase, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(["LB_WALBERLA"])
class LBBoundariesWalberlaSinglePrecision(LBBoundariesBase, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == "__main__":
    ut.main()
