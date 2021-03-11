#
# Copyright (C) 2010-2019 The ESPResSo project
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
import numpy as np

import espressomd.shapes


class ShapeTests(ut.TestCase):
    def test_Union(self):
        union = espressomd.shapes.Union()
        wall1 = espressomd.shapes.Wall(normal=[0, 0, 1], dist=0)
        wall2 = espressomd.shapes.Wall(normal=[0, 0, -1], dist=-10)
        self.assertTrue(union.call_method('empty'))
        union.add([wall1, wall2])
        self.assertFalse(union.call_method('empty'))
        self.assertEqual(union.size(), 2)

        # check object retrieval
        pwall1, pwall2 = union.call_method('get_elements')
        self.assertIsInstance(pwall1, espressomd.shapes.Wall)
        self.assertIsInstance(pwall2, espressomd.shapes.Wall)
        np.testing.assert_almost_equal(
            np.copy(pwall1.normal), np.copy(wall1.normal))
        np.testing.assert_almost_equal(
            np.copy(pwall2.normal), np.copy(wall2.normal))
        np.testing.assert_almost_equal(pwall1.dist, wall1.dist)
        np.testing.assert_almost_equal(pwall2.dist, wall2.dist)

        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 4.5])[0], 4.5)
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 5.0])[0], 5.0)
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 6.5])[0], 3.5)

        # negative distances are not well-defined for a union of shapes
        with self.assertRaises(ValueError):
            union.calc_distance(position=[1, 2, 11.5])
        union.clear()
        self.assertTrue(union.call_method('empty'))
        self.assertEqual(union.size(), 0)
        self.assertEqual(union.calc_distance(position=[1, 2, 6.5])[0], np.inf)

        union.add([wall1, wall2])
        union.remove(wall2)
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 6.5])[0], 6.5)


if __name__ == "__main__":
    ut.main()
