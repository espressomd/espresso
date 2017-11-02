#
# Copyright (C) 2017 The ESPResSo project
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

@ut.skipIf(not espressomd.has_features("LENNARD_JONES"),
           "Skipped because of LENNARD_JONES ")
class ForceCap(ut.TestCase):
    """Tests the force capping mechanism.

    """

    s = espressomd.System()
    s.cell_system.skin = 0.0
    s.seed = range(s.cell_system.get_state()["n_nodes"])

    @classmethod
    def setUpClass(cls):
        np.random.seed(42)

    def calc_f_max(self):
        f = np.power(self.s.part[:].f, 2)
        sqr_sum = (np.sum(f, axis=1))
        f_max = np.max(sqr_sum)**0.5
        return f_max

    def test(self):
        N = 200
        f_cap = 10.
        s = self.s
        s.part.clear()
        s.time_step = 0.1
        s.box_l = 3*[10.]
        s.part.add(pos=10. * np.random.random((N, 3)))

        self.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1000., sigma=2.,
            cutoff=1.5, shift=0.0)

        self.s.integrator.run(0)
        # Check that there is sth to cap
        self.assertGreater(self.calc_f_max(), f_cap)

        self.s.force_cap = f_cap

        # Check interface
        self.assertEqual(self.s.force_cap, f_cap)
        self.s.integrator.run(0)

        # Since there was a force larger than f_cap, the
        # maximum should now be f_cap.
        self.assertAlmostEqual(self.calc_f_max(), f_cap, places=7)

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
