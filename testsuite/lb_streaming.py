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
VISC = 20000000
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


class LBStreamingCommon(object):
    """
    Check the streaming step of the LB fluid implementation by setting all populations
    to zero except one. Relaxation is supressed by choosing appropriate parameters.

    """
    lbf = None
    system = espressomd.System(box_l=[3.0]*3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TAU
    grid = np.array([int(system.box_l[0] / AGRID), int(system.box_l[1] / AGRID), int(system.box_l[2] / AGRID)])

    def prepare(self):
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        self.reset_fluid_populations()

    def reset_fluid_populations(self):
        """Set all populations to 0.0.

        """
        for i in itertools.product(range(self.grid[0]), range(self.grid[1]), range(self.grid[2])):
            self.lbf[i].population = np.zeros(19)

    def set_fluid_populations(self, grid_index, n_v):
        """Set the population of direction n_v of grid_index to 1.0.

        """
        pop = np.zeros(19)
        pop[n_v] = 1.0
        self.lbf[grid_index].population = pop

    def test_population_streaming(self):
        self.prepare()
        for grid_index in itertools.product(range(self.grid[0]), range(self.grid[1]), range(self.grid[2])):
            for n_v in range(19):
                self.set_fluid_populations(grid_index,n_v)
                self.system.integrator.run(1)
                target_node = np.mod(grid_index + VELOCITY_VECTORS[n_v], self.grid[0])
                np.testing.assert_almost_equal(self.lbf[grid_index+VELOCITY_VECTORS[n_v]].population[n_v], 1.0)
                self.reset_fluid_populations()

@ut.skipIf(not espressomd.has_features('LB') or espressomd.has_features('SHANCHEN'),
        "Skipping test due to missing features.")
class LBCPU(ut.TestCase, LBStreamingCommon):
    """Test for the CPU implementation of the LB."""
    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(agrid=AGRID, tau=TAU, visc=VISC, bulk_visc=BULK_VISC, gamma_odd=1.0, gamma_even=1.0, fric=1.0, dens=1.0)

@ut.skipIf(not espressomd.has_features('LB_GPU') or espressomd.has_features('SHANCHEN'),
        "Skipping test due to missing features.")
class LBGPU(ut.TestCase, LBStreamingCommon):
    """Test for the GPU implementation of the LB."""
    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(agrid=AGRID, tau=TAU, visc=VISC, bulk_visc=BULK_VISC, gamma_odd=1.0, gamma_even=1.0, fric=1.0, dens=1.0)


if __name__ == "__main__":
    ut.main()
