#
# Copyright (C) 2010-2026 The ESPResSo project
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
BOUNDARY_VELOCITY = np.array([.0, 0.4, 0.5])
LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'kinematic_viscosity': VISC,
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
        self.system.lb = self.lbf

    def tearDown(self):
        self.system.lb = None

    def test(self):
        """
        Integrate the LB fluid until steady state is reached within a certain
        accuracy. Then compare the force balance between force exerted on fluid
        and forces acting on the boundaries.

        """
        wall_shape1 = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - AGRID))

        self.lbf.add_boundary_from_shape(wall_shape1, BOUNDARY_VELOCITY)
        self.lbf.add_boundary_from_shape(wall_shape2, BOUNDARY_VELOCITY)
        fluid_nodes = np.sum(np.logical_not(
            self.lbf[:, :, :].is_boundary).astype(int))

        self.system.integrator.run(20)
        diff = float("inf")
        old_val = float("inf")
        while diff > 0.00002:
            self.system.integrator.run(10)
            new_val = self.lbf.get_boundary_force_from_shape(wall_shape1)[0]
            diff = abs(new_val - old_val)
            old_val = new_val

        expected_force = np.copy(fluid_nodes * AGRID**3 * EXT_FORCE)

        measured_force_all = np.array(self.lbf.boundary_force)
        measured_force_1 = np.array(
            self.lbf.get_boundary_force_from_shape(wall_shape1))
        measured_force_2 = np.array(
            self.lbf.get_boundary_force_from_shape(wall_shape2))

        np.testing.assert_allclose(
            measured_force_all,
            expected_force,
            rtol=2E-2)

        np.testing.assert_allclose(
            measured_force_1,
            expected_force / 2,
            rtol=2E-2)

        np.testing.assert_allclose(
            measured_force_2,
            expected_force / 2,
            rtol=2E-2)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBForceWalberlaDoublePrecisionCPU(LBBoundaryForceCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBForceWalberlaSinglePrecisionCPU(LBBoundaryForceCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBBForceWalberlaDoublePrecisionGPU(LBBoundaryForceCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBBForceWalberlaSinglePrecisionGPU(LBBoundaryForceCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": True}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBForceWalberlaDoublePrecisionBlocks(
        LBBoundaryForceCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False,
                 "blocks_per_mpi_rank": [2, 2, 2]}


if __name__ == '__main__':
    ut.main()
