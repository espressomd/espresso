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
        with self.assertRaises(ValueError):
            c.shape = Wall(thisparameterdoesnotexist=0)

if __name__ == "__main__":
    ut.main()
