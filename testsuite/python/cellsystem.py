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
import tests_common


class CellSystem(ut.TestCase):
    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    n_nodes = system.cell_system.get_state()['n_nodes']

    def test_cell_system(self):
        parameters = {
            "n_square": {"use_verlet_lists": False},
            "regular_decomposition": {"use_verlet_lists": True},
            "hybrid_decomposition": {"use_verlet_lists": False,
                                     "n_square_types": {1, 3, 5},
                                     "cutoff_regular": 1.27},
        }
        for cell_system, params_in in parameters.items():
            setter = getattr(self.system.cell_system, f"set_{cell_system}")
            setter(**params_in)
            params_in["type"] = cell_system
            params_out = self.system.cell_system.get_state()
            tests_common.assert_params_match(self, params_in, params_out)
            params_out = self.system.cell_system.__getstate__()
            tests_common.assert_params_match(self, params_in, params_out)

    @ut.skipIf(n_nodes == 1, "Skipping test: only runs for n_nodes >= 2")
    def check_node_grid(self):
        for i in range(3):
            node_grid_ref = [1, 1, 1]
            node_grid_ref[i] = self.n_nodes
            self.system.cell_system.node_grid = node_grid_ref
            node_grid = self.system.cell_system.get_state()['node_grid']
            np.testing.assert_array_equal(node_grid, node_grid_ref)

    def test_node_grid_regular(self):
        self.system.cell_system.set_regular_decomposition()
        self.check_node_grid()

    def test_node_grid_hybrid(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=0)
        self.check_node_grid()


if __name__ == "__main__":
    ut.main()
