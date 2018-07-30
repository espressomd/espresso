
#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np


class Layered(ut.TestCase):
    S = espressomd.System(box_l=[1.0, 1.0, 1.0])
    S.seed = S.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.S.part.clear()
        self.S.cell_system.set_layered()

    def test_resort(self):
        n_part = 2351

        # Add the particles on node 0, so that they have to be
        # resorted
        for i in range(n_part):
            self.S.part.add(id=i, pos=[0, 0, 0], type=1)

        # And now change their positions
        for i in range(n_part):
            self.S.part[i].pos = pos = np.random.random(3)

        # Distribute the particles on the nodes
        part_dist = self.S.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # is still in a valid state after the particle exchange
        self.assertEqual(sum(self.S.part[:].type), n_part)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
