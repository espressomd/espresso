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
    'viscosity': VISC,
    'tau': TAU,
    'density': 1.0,
}


class LBStreamingCommon:

    """
    Check the streaming and relaxation steps of the LB fluid implementation by
    setting all populations to zero except one.

    """
    system = espressomd.System(box_l=[3.0] * 3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TAU
    grid = np.array(system.box_l / AGRID, dtype=int)

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMETERS, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()

    def reset_fluid_populations(self):
        """Set all populations to 0.0.

        """
        for i in itertools.product(range(self.grid[0]), range(
                self.grid[1]), range(self.grid[2])):
            self.lbf[i].population = np.zeros(19) + 1e-10

    def set_fluid_populations(self, grid_index):
        """Set the population of direction n_v of grid_index to n_v+1.

        """
        pop = np.arange(1, 20, dtype=float)
        self.lbf[grid_index].population = pop

    def test_population_streaming(self):
        self.reset_fluid_populations()
        for grid_index in itertools.product(
                range(self.grid[0]), range(self.grid[1]), range(self.grid[2])):
            self.set_fluid_populations(grid_index)
            self.system.integrator.run(1)
            for n_v in range(19):
                target_node_index = np.mod(
                    grid_index + VELOCITY_VECTORS[n_v], self.grid)
                np.testing.assert_allclose(
                    self.lbf[target_node_index].population[n_v],
                    REFERENCE_POPULATIONS[n_v], rtol=self.rtol)
                self.lbf[target_node_index].population = np.zeros(19) + 1e-10


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBStreamingWalberla(LBStreamingCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    rtol = 1e-10


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBStreamingWalberlaSinglePrecision(LBStreamingCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    rtol = 1e-5


# TODO WALBERLA
# @utx.skipIfMissingGPU()
# @utx.skipIfMissingFeatures(["WALBERLA"])
# class LBGPU(LBStreamingCommon, ut.TestCase):

#    """Test for the Walberla implementation of the LB on the GPU."""

#    lb_class = espressomd.lb.LBFluidWalberlaGPU
#    lb_params = {}
#    rtol = 1e-7


if __name__ == "__main__":
    ut.main()
