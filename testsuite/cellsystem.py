
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


class CellSystem(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def test_cell_system(self):
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        s = self.system.cell_system.get_state()
        self.assertEqual([s['use_verlet_list'], s['type']], [0, "nsquare"])
        self.system.cell_system.set_layered(n_layers=5)
        s = self.system.cell_system.get_state()
        self.assertEqual([s['type'], s['n_layers']], ["layered", 5])
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        s = self.system.cell_system.get_state()
        self.assertEqual(
            [s['use_verlet_list'], s['type']], [1, "domain_decomposition"])


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
