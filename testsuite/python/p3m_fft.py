#
# Copyright (C) 2020-2021 The ESPResSo project
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

import espressomd
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common

P3M_PARAMS = [
    {'cao': 7, 'r_cut': 3.103065490722656, 'alpha': 1.228153768561588, 'mesh': 48},
    {'cao': 7, 'r_cut': 4.477272033691406, 'alpha': 0.845808585620971, 'mesh': 32},
    {'cao': 7, 'r_cut': 2.393871545791626, 'alpha': 1.599093835130641, 'mesh': 64},
]

FFT_PLANS = {
    1: [([1, 1, 1], P3M_PARAMS[1])],
    2: [([2, 1, 1], P3M_PARAMS[1])],
    3: [([3, 1, 1], P3M_PARAMS[0])],
    4: [([2, 2, 1], P3M_PARAMS[1]),
        ([4, 1, 1], P3M_PARAMS[2])],
    6: [([3, 2, 1], P3M_PARAMS[0])],
    8: [([2, 2, 2], P3M_PARAMS[1]),
        ([4, 2, 1], P3M_PARAMS[2])],
}


@utx.skipIfMissingFeatures(["LENNARD_JONES", "P3M"])
class FFT_test(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    original_node_grid = tuple(system.cell_system.node_grid)
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.cell_system.node_grid = self.original_node_grid
        self.system.time_step = 0.01

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()

    def minimize(self):
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2**(1.0 / 6.0), shift="auto")
        self.system.integrator.set_steepest_descent(
            f_max=1, gamma=0.01, max_displacement=0.01)
        self.system.integrator.run(100)
        self.system.integrator.set_vv()
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0.0, sigma=1.0, cutoff=2)

    def add_charged_particles(self):
        np.random.seed(seed=42)
        num_pairs = 200
        positions = np.random.random((2 * num_pairs, 3))
        self.system.part.add(pos=positions * self.system.box_l,
                             q=num_pairs * [-1, 1])
        self.minimize()

    def add_magnetic_particles(self):
        np.random.seed(seed=42)
        num_part = 200
        positions = np.random.random((num_part, 3))
        dipoles = tests_common.random_dipoles(num_part)
        self.system.part.add(pos=positions * self.system.box_l,
                             dip=dipoles, rotation=num_part * [(1, 1, 1)])
        self.minimize()

    @ut.skipIf(n_nodes not in FFT_PLANS, f"no FFT plan for {n_nodes} threads")
    def test_fft_plans(self):
        import espressomd.electrostatics
        self.system.time_step = 0.01
        self.add_charged_particles()
        for node_grid, p3m_params in FFT_PLANS[self.n_nodes]:
            self.system.cell_system.node_grid = node_grid
            solver = espressomd.electrostatics.P3M(
                prefactor=2, accuracy=1e-6, tune=False, **p3m_params)
            self.system.actors.add(solver)
            ref_energy = -75.871906
            p3m_energy = self.system.analysis.energy()['coulomb']
            self.system.actors.clear()
            np.testing.assert_allclose(p3m_energy, ref_energy, rtol=1e-4)

    @utx.skipIfMissingFeatures("P3M")
    @ut.skipIf(n_nodes < 2 or n_nodes >= 8, "only runs for 2 <= n_nodes <= 7")
    def test_unsorted_node_grid_exception_p3m(self):
        import espressomd.electrostatics
        self.system.time_step = 0.01
        self.add_charged_particles()
        unsorted_node_grid = self.system.cell_system.node_grid[::-1]
        self.system.cell_system.node_grid = unsorted_node_grid
        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: P3M_init: node grid must be sorted, largest first'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    @ut.skipIf(n_nodes < 2 or n_nodes >= 8, "only runs for 2 <= n_nodes <= 7")
    def test_unsorted_node_grid_exception_dp3m(self):
        import espressomd.magnetostatics
        self.system.time_step = 0.01
        self.add_magnetic_particles()
        unsorted_node_grid = self.system.cell_system.node_grid[::-1]
        self.system.cell_system.node_grid = unsorted_node_grid
        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: dipolar P3M_init: node grid must be sorted, largest first'):
            self.system.actors.add(solver)


if __name__ == "__main__":
    ut.main()
