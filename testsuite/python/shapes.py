import unittest as ut
import numpy as np

import espressomd.shapes


class ShapeTests(ut.TestCase):
    def test_Union(self):
        union = espressomd.shapes.Union()
        wall1 = espressomd.shapes.Wall(normal=[0, 0, 1], dist=0)
        wall2 = espressomd.shapes.Wall(normal=[0, 0, -1], dist=-10)
        union.add([wall1, wall2])
        self.assertEqual(union.size(), 2)

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
        self.assertEqual(union.size(), 0)
        self.assertEqual(union.calc_distance(position=[1, 2, 6.5])[0], np.inf)

        union.add([wall1, wall2])
        union.remove(wall2)
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 6.5])[0], 6.5)


if __name__ == "__main__":
    ut.main()
