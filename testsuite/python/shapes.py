#
# Copyright (C) 2010-2026 The ESPResSo project
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

    def test_torus(self):
        """
        Check the torus shape. Nomenclature: phi for the toroidal direction
        and theta for the poloidal direction.
        """
        # check points on the revolution axis, using the intercept theorem
        shape = espressomd.shapes.Torus(
            normal=[0, 0, 1], direction=1, radius=1.5, tube_radius=0.1)
        for z in np.arange(-5, 6):
            dist, vec = shape.calc_distance(
                position=[0., 0., z * shape.radius])
            h = np.sqrt((z * shape.radius)**2 + shape.radius**2)
            np.testing.assert_allclose(dist, h - shape.tube_radius)
            np.testing.assert_allclose(
                np.copy(vec), dist / h * shape.radius * np.array([-1., 0., z]))

        shape = espressomd.shapes.Torus(center=[15] * 3, normal=[0, 0, 1],
                                        direction=1, radius=10., tube_radius=1.)
        # check points on the tube surface, at different altitudes
        for phi in np.linspace(0., 2. * np.pi, 65):
            origin = np.copy(shape.center) + [shape.radius, 0., 0.]
            for r in np.linspace(0., 2. * shape.tube_radius, 11)[1:]:
                if np.abs(r - shape.tube_radius) < 1e-6:
                    continue
                ref_dist = shape.tube_radius - r
                for theta in np.linspace(0., 2. * np.pi, 65):
                    unit_vec = np.array([np.cos(theta), 0., np.sin(theta)])
                    dist, vec = shape.calc_distance(
                        position=origin + r * unit_vec)
                    np.testing.assert_allclose(dist, -ref_dist)
                    np.testing.assert_allclose(
                        np.copy(vec), -ref_dist * unit_vec, atol=1e-10)
            # check points exactly on the tube surface
            r = shape.tube_radius
            for theta in np.linspace(0., 2. * np.pi, 65):
                unit_vec = np.array([np.cos(theta), 0., np.sin(theta)])
                dist, vec = shape.calc_distance(position=origin + r * unit_vec)
                np.testing.assert_allclose(dist, 0., atol=1e-14)
                np.testing.assert_allclose(np.copy(vec), 0., atol=1e-14)

        # check points on the core circle; angle is undefined and by default
        # the vector will point toward the torus center
        for phi in np.linspace(0., 2. * np.pi, 65):
            unit_vec = np.array([np.cos(phi), np.sin(phi), 0.])
            dist, vec = shape.calc_distance(
                position=shape.center + shape.radius * unit_vec)
            ref_vec = -shape.tube_radius * unit_vec
            np.testing.assert_allclose(np.copy(vec), ref_vec, atol=1e-14)


if __name__ == "__main__":
    ut.main()
