#
# Copyright (C) 2017-2018 The ESPResSo project
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
import unittest_decorators as utx
import espressomd


@utx.skipIfMissingFeatures("EXCLUSIONS")
class AutoExclusions(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.part.clear()

    def test_linear(self):
        bond = espressomd.interactions.Virtual()
        s = self.system
        s.bonded_inter.add(bond)

        for i in range(10):
            s.part.add(id=i, pos=[0, 0, 0])

        for i in range(9):
            s.part[i].add_bond((bond, i + 1))

        s.auto_exclusions(1)

        for p in range(1, 9):
            excl = s.part[p].exclusions
            self.assertEqual(len(excl), 2)
            self.assertTrue((p - 1) in excl)
            self.assertTrue((p + 1) in excl)

        excl = s.part[0].exclusions
        self.assertEqual(len(excl), 1)
        self.assertTrue(1 in excl)

        excl = s.part[9].exclusions
        self.assertEqual(len(excl), 1)
        self.assertTrue(8 in excl)

    def test_ring(self):
        bond = espressomd.interactions.Virtual()
        s = self.system
        s.bonded_inter.add(bond)

        for i in range(10):
            s.part.add(id=i, pos=[0, 0, 0])

        for i in range(10):
            s.part[i].add_bond((bond, (i + 1) % 10))

        s.auto_exclusions(2)

        for p in range(10):
            excl = s.part[p].exclusions
            self.assertEqual(len(excl), 4)
            for i in range(1, 3):
                self.assertTrue((p - i) % 10 in excl)
                self.assertTrue((p + i) % 10 in excl)

if __name__ == "__main__":
    ut.main()
