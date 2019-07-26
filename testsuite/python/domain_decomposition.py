#
# Copyright (C) 2013-2018 The ESPResSo project
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


class DomainDecomposition(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    def setUp(self):
        self.system.part.clear()
        self.system.cell_system.set_domain_decomposition(
            use_verlet_lists=False)

    def test_resort(self):
        n_part = 2351

        # Add the particles on node 0, so that they have to be resorted
        for i in range(n_part):
            self.system.part.add(id=i, pos=[0, 0, 0], type=1)

        # And now change their positions
        for i in range(n_part):
            self.system.part[i].pos = pos = np.random.random(3)

        # Distribute the particles on the nodes
        part_dist = self.system.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # is still in a valid state after the particle exchange
        self.assertEqual(sum(self.system.part[:].type), n_part)

    def test_min_num_cells(self):
        s = self.system
        cs = s.cell_system
        cs.min_num_cells = 23

        self.assertEqual(cs.min_num_cells, 23)
        cell_grid = cs.get_state()['cell_grid']
        n_cells = cell_grid[0] * cell_grid[1] * cell_grid[2]
        # Check that we have neither too few nor too many cells
        self.assertGreaterEqual(n_cells, cs.min_num_cells)
        self.assertLessEqual(n_cells, cs.max_num_cells)

if __name__ == "__main__":
    ut.main()
