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
from __future__ import print_function
import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes
import unittest as ut
import numpy as np


@ut.skipIf(not espressomd.has_features(["LB", "LB_BOUNDARIES"]),
           "Features not available, skipping test.")
class LBBoundaryVelocityTest(ut.TestCase):

    """Test slip velocity of boundaries.

       In this simple test add wall with a slip verlocity is
       added and checkeckt if the fluid obtains the same velocity.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = .5
    system.cell_system.skin = 0.1

    def test(self):
        system = self.system

        lb_fluid = espressomd.lb.LBFluid(
            agrid=2.0, dens=1.0, visc=1.0, tau=0.03)
        system.actors.add(lb_fluid)

        v_boundary = [0.03, 0.02, 0.01]

        wall_shape = espressomd.shapes.Wall(normal=[1, 2, 3], dist=0.5)
        wall = espressomd.lbboundaries.LBBoundary(
            shape=wall_shape, velocity=v_boundary)
        system.lbboundaries.add(wall)

        system.integrator.run(10000)

        v_fluid = lb_fluid[1, 0, 0].velocity
        self.assertAlmostEqual(v_fluid[0], v_boundary[0], places=3)
        self.assertAlmostEqual(v_fluid[1], v_boundary[1], places=3)
        self.assertAlmostEqual(v_fluid[2], v_boundary[2], places=3)

if __name__ == "__main__":
    ut.main()
