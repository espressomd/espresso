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

import espressomd.lb
import espressomd.shapes

AGRID = 0.5
EXT_FORCE = np.array([-.01, 0.02, 0.03])
VISC = 3.5
DENS = 1.5
TIME_STEP = 0.05
LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'viscosity': VISC,
             'tau': TIME_STEP,
             'ext_force_density': EXT_FORCE}


class LBBoundaryForceCommon:

    """
    Checks force on lb boundaries for a fluid with a uniform volume force
    """

    system = espressomd.System(box_l=np.array([12.0, 4.0, 4.0]) * AGRID)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()

    def test(self):
        """
        Integrate the LB fluid until steady state is reached within a certain
        accuracy. Then compare the force balance between force exerted on fluid
        and forces acting on the boundaries.

        """
        wall_shape1 = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - AGRID))

        fluid_nodes = np.sum(np.logical_not(
            self.lbf[:, :, :].is_boundary).astype(int))
        self.lbf.add_boundary_from_shape(wall_shape1)
        self.lbf.add_boundary_from_shape(wall_shape2)

        # TODO WALBERLA: (#4381)
        self.skipTest("boundary forces not implemented at the moment")

        self.system.integrator.run(20)
        diff = float("inf")
        old_val = float("inf")
        while diff > 0.002:
            self.system.integrator.run(10)
            new_val = self.lbf.boundary['wall1'].get_force()[0]
            diff = abs(new_val - old_val)
            old_val = new_val

        expected_force = fluid_nodes * AGRID**3 * \
            np.copy(self.lbf.ext_force_density)
        measured_force = np.array(self.lbf.boundary['wall1'].get_force()) + \
            np.array(self.lbf.boundary['wall2'].get_force())
        # TODO WALBERLA: the force converges to 90% of the expected force
        np.testing.assert_allclose(
            measured_force,
            expected_force * 0.9,
            atol=1E-10)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBoundaryForceWalberla(LBBoundaryForceCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBoundaryForceWalberlaSinglePrecision(
        LBBoundaryForceCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == '__main__':
    ut.main()
