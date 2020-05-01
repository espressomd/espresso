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


class DomainDecomposition(ut.TestCase):
    system = espressomd.System(box_l=[50.0, 50.0, 50.0])

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
            self.system.part[i].pos = self.system.box_l * np.random.random(3)

        # Distribute the particles on the nodes
        part_dist = self.system.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # is still in a valid state after the particle exchange
        self.assertEqual(sum(self.system.part[:].type), n_part)

    def test_position_rounding(self):
        """This places a particle on the box boundary,
           with parameters that could cause problems with
           rounding."""
        self.system.box_l = [50.0, 50.0, 50.0]
        self.system.cell_system.skin = 0.4
        self.system.min_global_cut = 12.0 / 4.25
        self.system.part.add(pos=[25, 25, 0])
        self.assertEqual(1, len(self.system.part))


if __name__ == "__main__":
    ut.main()
