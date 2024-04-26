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
import itertools
import numpy as np

import espressomd
import espressomd.lb

"""
Tests for the streaming of populations of the LB algorithm.

"""

AGRID = 0.5
TAU = 0.1
VISC = 1e18
VELOCITY_VECTORS = np.array([
    [0, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [-1, 0, 0],
    [1, 0, 0],
    [0, 0, 1],
    [0, 0, -1],
    [-1, 1, 0],
    [1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
    [0, 1, 1],
    [0, -1, 1],
    [-1, 0, 1],
    [1, 0, 1],
    [0, 1, -1],
    [0, -1, -1],
    [-1, 0, -1],
    [1, 0, -1]])
# populations after streaming and relaxation using parameters omega_odd = 2
# and omega_bulk = omega_even = omega_shear = 0
REFERENCE_POPULATIONS = np.array([
    1,
    2 / 3,
    4 + 1 / 3,
    3 + 1 / 3,
    5 + 2 / 3,
    1 + 1 / 3,
    11 + 2 / 3,
    9,
    9 + 2 / 3,
    9 + 1 / 3,
    10,
    13,
    14 + 1 / 3,
    15 + 1 / 3,
    16,
    14 + 2 / 3,
    16,
    17,
    17 + 2 / 3])
LB_PARAMETERS = {
    'agrid': AGRID,
    'kinematic_viscosity': VISC,
    'tau': TAU,
    'density': 1.0,
}


class LBStreamingCommon:

    """
    Check the streaming and relaxation steps of the LB fluid implementation by
    setting all populations to zero except for one cell.

    """
    system = espressomd.System(box_l=[3., 2., 2.])
    system.cell_system.skin = 0.1 * AGRID
    system.time_step = TAU

    def setUp(self):
        self.system.box_l = self.box_l
        self.lbf = self.lb_class(**LB_PARAMETERS, **self.lb_params)
        self.system.lb = self.lbf

    def tearDown(self):
        self.system.lb = None

    def test_population_streaming(self):
        pop_default = np.zeros(19) + 1e-10
        pop_source = np.arange(1, 20, dtype=float)
        grid = np.array(self.system.box_l / AGRID, dtype=int)

        # reset fluid populations
        self.lbf[:, :, :].population = pop_default

        # check streaming
        for grid_index in itertools.product(
                range(grid[0]), range(grid[1]), range(grid[2])):
            self.lbf[grid_index].population = pop_source
            self.system.integrator.run(1)
            for n_v in range(19):
                dst_vec = np.array(VELOCITY_VECTORS[n_v])
                target_node_index = np.mod(grid_index + dst_vec, grid)
                np.testing.assert_allclose(
                    self.lbf[target_node_index].population[n_v],
                    REFERENCE_POPULATIONS[n_v], rtol=self.rtol,
                    err_msg=f"streaming is incorrect in direction {VELOCITY_VECTORS[n_v]} from cell at {grid_index}")
                self.lbf[target_node_index].population = pop_default


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBStreamingWalberlaDoublePrecisionCPU(LBStreamingCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    box_l = [3., 2., 2.]
    rtol = 1e-10


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBStreamingWalberlaSinglePrecisionCPU(LBStreamingCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    box_l = [3., 2., 2.]
    rtol = 1e-5


@utx.skipIfMissingGPU()
@ut.skipIf(LBStreamingCommon.system.cell_system.get_state()["n_nodes"] > 2,
           "only runs for 2 or less MPI ranks")
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBStreamingWalberlaDoublePrecisionGPU(LBStreamingCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": False}
    box_l = [2., 1.5, 1.5]
    rtol = 1e-10


@utx.skipIfMissingGPU()
@ut.skipIf(LBStreamingCommon.system.cell_system.get_state()["n_nodes"] > 2,
           "only runs for 2 or less MPI ranks")
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBStreamingWalberlaSinglePrecisionGPU(LBStreamingCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": True}
    box_l = [2., 1.5, 1.5]
    rtol = 1e-5


if __name__ == "__main__":
    ut.main()
