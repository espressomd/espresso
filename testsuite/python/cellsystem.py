
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
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond


class CellSystem(ut.TestCase):
    S = espressomd.System()

    def test_cell_system(self):
        self.S.cell_system.set_n_square(use_verlet_lists=False)
        s = self.S.cell_system.get_state()
        self.assertEqual(s, {"use_verlet_lists": 0, "type": "nsquare"})
        self.S.cell_system.set_layered(nLayers=5)
        s = self.S.cell_system.get_state()
        self.assertEqual(s, {"type": "layered", "nLayers": 5})

        self.S.cell_system.set_domain_decomposition(use_verlet_lists=True)
        s = self.S.cell_system.get_state()
        self.assertEqual(
            s, {"use_verlet_lists": 1, "type": "domain_decomposition"})


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
