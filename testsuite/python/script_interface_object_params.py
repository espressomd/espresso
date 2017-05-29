from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd.shapes import Wall,Sphere
from espressomd.constraints import ShapeBasedConstraint

@ut.skipIf(not espressomd.has_features(["CONSTRAINTS"]),"Tests requries CONSTRAINTS")
class ScriptInterfaceObjectParams(ut.TestCase):
    """Tests that object parameters are assigned the correct python class"""

    def test(self):
        c=ShapeBasedConstraint()
        w=Wall(normal=[-1,0,0])
        c.shape=w
        # Does the shape parameter return the correct lcass
        self.assertEqual(c.shape.__class__,Wall)
        # Do the sciprt object match
        self.assertTrue(c.shape==w)
        
        # Different shape
        c.shape=Sphere(radius=1)
        # Test class
        self.assertTrue(c.shape.__class__==Sphere)
        # Test parameter retrieval
        self.assertTrue(abs(c.shape.radius-1)<=1E-8)




if __name__ == "__main__":
    ut.main()
