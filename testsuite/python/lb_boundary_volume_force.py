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

import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes

from tests_common import count_fluid_nodes


"""
Checks force on lb boundaries for a fluid with a uniform volume force
"""


AGRID = 0.5
EXT_FORCE = np.array([-.01, 0.02, 0.03])
VISC = 3.5
DENS = 1.5
TIME_STEP = 0.05
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             'ext_force_density': EXT_FORCE}


class LBBoundaryForceCommon:

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=np.array([12.0, 4.0, 4.0]) * AGRID)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def test(self):
        """
        Integrate the LB fluid until steady state is reached within a certain
        accuracy. Then compare the force balance between force exerted on fluid
        and forces acting on the boundaries.

        """
        self.system.actors.clear()
        self.system.lbboundaries.clear()
        self.system.actors.add(self.lbf)
        wall_shape1 = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - AGRID))
        wall1 = espressomd.lbboundaries.LBBoundary(shape=wall_shape1)
        wall2 = espressomd.lbboundaries.LBBoundary(shape=wall_shape2)

        self.system.lbboundaries.add(wall1)
        self.system.lbboundaries.add(wall2)
        fluid_nodes = count_fluid_nodes(self.lbf)

        self.system.integrator.run(20)
        diff = float("inf")
        old_val = float("inf")
        while diff > 0.002:
            self.system.integrator.run(10)
            new_val = wall1.get_force()[0]
            diff = abs(new_val - old_val)
            old_val = new_val

        surface_area = self.system.box_l[1] * self.system.box_l[2]
        expected_force = fluid_nodes * AGRID**3 * \
            np.copy(self.lbf.ext_force_density)
        measured_force = np.array(wall1.get_force()) + \
            np.array(wall2.get_force())
        np.testing.assert_allclose(measured_force, expected_force, atol=2E-2)


@utx.skipIfMissingFeatures(['LB_BOUNDARIES', 'EXTERNAL_FORCES'])
class LBCPUBoundaryForce(ut.TestCase, LBBoundaryForceCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(['LB_BOUNDARIES_GPU', 'EXTERNAL_FORCES'])
class LBGPUBoundaryForce(ut.TestCase, LBBoundaryForceCommon):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
