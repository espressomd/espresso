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

    def test_cell_system(self):
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        s = self.system.cell_system.get_state()
        self.assertEqual([s['use_verlet_list'], s['type']], [0, "nsquare"])
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        s = self.system.cell_system.get_state()
        self.assertEqual(
            [s['use_verlet_list'], s['type']], [1, "domain_decomposition"])

    def test_node_grid(self):
        self.system.cell_system.set_domain_decomposition()
        n_nodes = self.system.cell_system.get_state()['n_nodes']
        if n_nodes == 1:
            return
        self.system.cell_system.node_grid = [n_nodes, 1, 1]
        s = self.system.cell_system.get_state()
        np.testing.assert_array_equal(
            s['node_grid'], [n_nodes, 1, 1])

    def test_fully_connected_node_grid(self):
        n_nodes = self.system.cell_system.get_state()['n_nodes']
        if n_nodes == 1:
            return
        self.system.cell_system.node_grid = [1, 1, n_nodes]
        with self.assertRaises(Exception):
            self.system.cell_system.set_domain_decomposition(
                fully_connected=[False, False, True])

    def test_particle_pair(self):
        n_nodes = self.system.cell_system.get_state()['n_nodes']
        if n_nodes == 1:
            return

        self.system.cell_system.node_grid = [1, 1, n_nodes]
        self.system.cell_system.set_domain_decomposition(
            fully_connected=[True, True, False])

        self.system.part.add(id=0, pos=[2.5, 4.75, 2.5])
        self.system.part.add(id=1, pos=[1.0, 5.25, 2.5])

        pairs = self.system.cell_system.get_pairs_(2.5)
        np.testing.assert_array_equal(pairs, [[0, 1]])


if __name__ == "__main__":
    ut.main()
