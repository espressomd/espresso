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
import unittest as ut

import espressomd.shapes
import espressomd.constraints


class ScriptInterfaceObjectParams(ut.TestCase):

    """Tests that object parameters are assigned the correct python class"""

    def test(self):
        c = espressomd.constraints.ShapeBasedConstraint()
        w = espressomd.shapes.Wall(normal=[-1, 0, 0])
        c.shape = w
        # Does the shape parameter return the correct class
        self.assertEqual(c.shape.__class__, espressomd.shapes.Wall)
        # Does the script object match
        self.assertEqual(c.shape, w)

        # Different shape
        c.shape = espressomd.shapes.Sphere(radius=1)
        # Test class
        self.assertEqual(c.shape.__class__, espressomd.shapes.Sphere)
        # Test parameter retrieval
        self.assertAlmostEqual(c.shape.radius, 1, places=8)
        with self.assertRaises(RuntimeError):
            c.shape = espressomd.shapes.Wall(thisparameterdoesnotexist=0)

    def test_compare(self):
        a = espressomd.constraints.ShapeBasedConstraint()
        b = espressomd.constraints.ShapeBasedConstraint()
        c = 5
        self.assertEqual(a, a)
        self.assertNotEqual(a, b)
        self.assertFalse(a == c)
        self.assertFalse(c == a)
        self.assertTrue(a != c)
        self.assertTrue(c != a)
        with self.assertRaises(NotImplementedError):
            a > a
        with self.assertRaises(NotImplementedError):
            a >= a
        with self.assertRaises(NotImplementedError):
            a < a
        with self.assertRaises(NotImplementedError):
            a <= a
        with self.assertRaises(NotImplementedError):
            a > c
        with self.assertRaises(NotImplementedError):
            a >= c
        with self.assertRaises(NotImplementedError):
            a < c
        with self.assertRaises(NotImplementedError):
            a <= c


if __name__ == "__main__":
    ut.main()
