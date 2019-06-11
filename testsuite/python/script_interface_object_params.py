# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd.shapes import Wall, Sphere
from espressomd.constraints import ShapeBasedConstraint


class ScriptInterfaceObjectParams(ut.TestCase):

    """Tests that object parameters are assigned the correct python class"""

    def test(self):
        c = ShapeBasedConstraint()
        w = Wall(normal=[-1, 0, 0])
        c.shape = w
        # Does the shape parameter return the correct lcass
        self.assertEqual(c.shape.__class__, Wall)
        # Do the sciprt object match
        self.assertEqual(c.shape, w)

        # Different shape
        c.shape = Sphere(radius=1)
        # Test class
        self.assertEqual(c.shape.__class__, Sphere)
        # Test parameter retrieval
        self.assertAlmostEqual(c.shape.radius, 1, places=8)
        with self.assertRaises(RuntimeError):
            c.shape = Wall(thisparameterdoesnotexist=0)

if __name__ == "__main__":
    ut.main()
