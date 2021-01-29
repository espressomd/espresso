#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd
import numpy as np


class CellSystem(ut.TestCase):
    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    n_nodes = system.cell_system.get_state()['n_nodes']

    def test_cell_system(self):
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        s = self.system.cell_system.get_state()
        self.assertEqual([s['use_verlet_list'], s['type']], [0, "nsquare"])
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        s = self.system.cell_system.get_state()
        self.assertEqual(
            [s['use_verlet_list'], s['type']], [1, "domain_decomposition"])

    @ut.skipIf(n_nodes == 1, "Skipping test: only runs for n_nodes >= 2")
    def test_node_grid(self):
        self.system.cell_system.set_domain_decomposition()
        for i in range(3):
            node_grid_ref = [1, 1, 1]
            node_grid_ref[i] = self.n_nodes
            self.system.cell_system.node_grid = node_grid_ref
            node_grid = self.system.cell_system.get_state()['node_grid']
            np.testing.assert_array_equal(node_grid, node_grid_ref)


if __name__ == "__main__":
    ut.main()
