#
# Copyright (C) 2013-2022 The ESPResSo project
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
from espressomd.interactions import FeneBond
import numpy as np
import unittest as ut


class CutOff(ut.TestCase):

    """
    Test interaction cutoffs
    """

    def test(self):
        system = espressomd.System(box_l=[15.0, 15.0, 15.0])
        system.cell_system.skin = 1

        # Initial state. Skin does not influence cutoffs as long as there are
        # no interactions
        self.assertEqual(system.cell_system.max_cut_nonbonded, -1)
        self.assertEqual(system.cell_system.max_cut_bonded, -1)
        self.assertEqual(system.cell_system.interaction_range, -1)

        # Bonded interaction
        fene = FeneBond(r_0=1, d_r_max=2, k=1)
        system.bonded_inter.add(fene)
        self.assertEqual(system.cell_system.max_cut_bonded, 3)
        n_nodes = np.prod(system.cell_system.node_grid)
        if n_nodes == 1:
            # Bonds don't influence interaction range
            self.assertEqual(system.cell_system.interaction_range, -1)
        else:
            self.assertEqual(system.cell_system.interaction_range,
                             system.cell_system.max_cut_bonded +
                             system.cell_system.skin)

        system.bonded_inter.remove(fene._bond_id)
        self.assertEqual(system.cell_system.max_cut_bonded, -1)
        self.assertEqual(system.cell_system.interaction_range, -1)

        if espressomd.has_features("LENNARD_JONES"):
            lj_off_params = system.non_bonded_inter[0,
                                                    0].lennard_jones.get_params()
            system.non_bonded_inter[0, 0].lennard_jones.set_params(
                sigma=1, epsilon=1, cutoff=2.5, shift="auto")
            self.assertEqual(
                system.cell_system.max_cut_nonbonded, 2.5)
            self.assertEqual(system.cell_system.interaction_range,
                             system.cell_system.max_cut_nonbonded +
                             system.cell_system.skin)

            system.non_bonded_inter[0,
                                    0].lennard_jones.set_params(**lj_off_params)
            self.assertEqual(system.cell_system.max_cut_nonbonded, -1)
            self.assertEqual(system.cell_system.interaction_range, -1)


if __name__ == '__main__':
    ut.main()
