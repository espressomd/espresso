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


@ut.skipIf(not espressomd.has_features(["VIRTUAL_SITES"]),
           "Features not available, skipping test.")
class LBBoundaryThermoVirtualTest(ut.TestCase):

    """Test slip velocity of boundaries.

       In this simple test add wall with a slip verlocity is
       added and checkeckt if the fluid obtains the same velocity.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    system.time_step = 1.0
    system.cell_system.skin = 0.1

    def tearDown(self):
        self.system.part.clear()

        for a in self.system.actors:
            self.system.actors.remove(a)

    def check_virtual(self, fluid_class):
        s = self.system
        lb_fluid = fluid_class(
            agrid=1.0, dens=1.0, visc=1.0, tau=1.0, kT=0.0)
        s.actors.add(lb_fluid)

        virtual = s.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = s.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        s.thermostat.set_lb(
            LB_fluid=lb_fluid,
            act_on_virtual=False,
            gamma=1.0)

        s.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])

        s.thermostat.set_lb(LB_fluid=lb_fluid, act_on_virtual=True)

        virtual.v = [1, 0, 0]
        physical.v = [1, 0, 0]

        s.actors.remove(lb_fluid)
        lb_fluid = fluid_class(
            agrid=1.0, dens=1.0, visc=1.0, tau=1.0)
        s.actors.add(lb_fluid)
        s.thermostat.set_lb(LB_fluid=lb_fluid, gamma=1.0)
        virtual.pos = physical.pos
        virtual.v = 1, 0, 0
        physical.v = 1, 0, 0
        s.integrator.run(1)

        # The forces are not exactly -1. because the fluid is not at
        # rest anymore because of the previous check.
        np.testing.assert_almost_equal(np.copy(physical.f), np.copy(virtual.f))
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])
        np.testing.assert_almost_equal(np.copy(virtual.f), [-1, 0, 0])

    def test_lb_cpu(self):
        self.check_virtual(espressomd.lb.LBFluid)

    @ut.skipIf(
        not espressomd.gpu_available() or not espressomd.has_features(
            ["CUDA"]),
               "Features or gpu not available, skipping test.")
    def test_lb_gpu(self):
        self.check_virtual(espressomd.lb.LBFluidGPU)

if __name__ == "__main__":
    ut.main()
