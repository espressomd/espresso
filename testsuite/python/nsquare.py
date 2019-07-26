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


class NSquare(ut.TestCase):
    S = espressomd.System(box_l=[1.0, 1.0, 1.0])
    S.seed = S.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.S.part.clear()
        self.S.cell_system.set_n_square(use_verlet_lists=False)

    def test_load_balancing(self):
        n_part = 235
        n_nodes = self.S.cell_system.get_state()['n_nodes']
        n_part_avg = n_part // n_nodes

        for i in range(n_part):
            self.S.part.add(id=i, pos=np.random.random(3), type=1)

        part_dist = self.S.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part)

        # Check that the particules are evenly distributed
        for node_parts in part_dist:
            self.assertLess(abs(node_parts - n_part_avg), 2)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # are still in a valid state after the particle exchange
        self.assertEqual(sum(self.S.part[:].type), n_part)


if __name__ == "__main__":
    ut.main()
