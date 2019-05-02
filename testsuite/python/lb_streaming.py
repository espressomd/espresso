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
BULK_VISC = VISC
VELOCITY_VECTORS = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, 1],
    [0, -1, -1],
    [0, 1, -1],
    [0, -1, 1]])
LB_PARAMETERS = {
    'agrid': AGRID,
    'visc': VISC,
    'bulk_visc': BULK_VISC,
    'tau': TAU,
    'dens': 1.0,
    'gamma_odd': 1.0,
    'gamma_even': 1.0
}


class LBStreamingCommon(object):

    """
    Check the streaming step of the LB fluid implementation by setting all populations
    to zero except one. Relaxation is supressed by choosing appropriate parameters.

    """
    lbf = None
    system = espressomd.System(box_l=[3.0] * 3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TAU
    grid = np.array([int(system.box_l[0] / AGRID),
                     int(system.box_l[1] / AGRID), int(system.box_l[2] / AGRID)])

    def prepare(self):
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        self.reset_fluid_populations()

    def reset_fluid_populations(self):
        """Set all populations to 0.0.

        """
        for i in itertools.product(range(self.grid[0]), range(self.grid[1]), range(self.grid[2])):
            self.lbf[i].population = np.zeros(19)

    def set_fluid_populations(self, grid_index):
        """Set the population of direction n_v of grid_index to n_v+1.

        """
        pop = np.arange(1, 20)
        self.lbf[grid_index].population = pop

    def test_population_streaming(self):
        self.prepare()
        for grid_index in itertools.product(range(self.grid[0]), range(self.grid[1]), range(self.grid[2])):
            self.set_fluid_populations(grid_index)
            self.system.integrator.run(1)
            for n_v in range(19):
                target_node_index = np.mod(
                    grid_index + VELOCITY_VECTORS[n_v], self.grid)
                np.testing.assert_almost_equal(
                    self.lbf[target_node_index].population[n_v], float(n_v + 1))
                self.lbf[target_node_index].population = np.zeros(19)


class LBCPU(ut.TestCase, LBStreamingCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMETERS)


@ut.skipIf(not espressomd.gpu_available() or 
           not espressomd.has_features(
    'LB_GPU'),
    "Skipping test due to missing features.")
class LBGPU(ut.TestCase, LBStreamingCommon):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMETERS)


if __name__ == "__main__":
    ut.main()
