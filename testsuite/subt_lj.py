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

from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np

@ut.skipIf(not espressomd.has_features("LENNARD_JONES"),
           "Skipped because of not LENNARD_JONES")
class SubtLjTest(ut.TestCase):
    system = espressomd.System(box_l=[10, 10, 10])
    system.time_step = .1

    def setUp(self):
        self.system.part.clear()

    def test(self):
        s = self.system
        s.part.add(id=0, pos=[4.5, 4.5, 4.5], type=0)
        s.part.add(id=1, pos=[5.5, 5.5, 5.5], type=0)

        s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=2.14, sigma=2.56,
            cutoff=1.5, shift=1.1, offset=0.5)

        s.integrator.run(0)
        f = np.sum(s.part[:].f**2)

        self.assertGreater(f, 10.)

        subt = espressomd.interactions.SubtLJ()
        s.bonded_inter.add(subt)

        s.part[0].add_bond((subt, 1))

        s.integrator.run(0)
        f = np.sum(s.part[:].f**2)

        print(s.analysis.energy())

        self.assertAlmostEqual(f, 0, places=10)
        self.assertAlmostEqual(s.analysis.energy()['total'] , 0, places=10)

if __name__ == "__main__":
    ut.main()
