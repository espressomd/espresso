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
import unittest_decorators as utx
import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.lbboundaries
from itertools import product


class LBBoundariesBase:
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.1

    wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
    wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-7.5)

    def test_add(self):
        boundary = espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1)

        self.system.lbboundaries.add(boundary)
        self.assertEqual(boundary, self.system.lbboundaries[0])

    def test_remove(self):
        lbb = self.system.lbboundaries

        b1 = lbb.add(
            espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        b2 = lbb.add(
            espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))

        lbb.remove(b1)

        self.assertNotIn(b1, lbb)
        self.assertIn(b2, lbb)

    def test_size(self):
        lbb = self.system.lbboundaries
        self.assertEqual(lbb.size(), 0)

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        self.assertEqual(lbb.size(), 1)

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        self.assertEqual(lbb.size(), 2)

    def test_empty(self):
        lbb = self.system.lbboundaries
        self.assertTrue(lbb.empty())

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        self.assertFalse(lbb.empty())

    def test_clear(self):
        lbb = self.system.lbboundaries

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))

        lbb.clear()

        self.assertTrue(lbb.empty())

    def check_boundary_flags(self, boundarynumbers):
        rng = range(20)

        for i in product(range(0, 5), rng, rng):
            self.assertEqual(self.lbf[i].boundary, boundarynumbers[0])

        for i in product(range(5, 15), rng, rng):
            self.assertEqual(self.lbf[i].boundary, boundarynumbers[1])

        for i in product(range(15, 20), rng, rng):
            self.assertEqual(self.lbf[i].boundary, boundarynumbers[2])

        self.system.lbboundaries.clear()
        for i in product(rng, rng, rng):
            self.assertEqual(self.lbf[i].boundary, 0)

    def test_boundary_flags(self):
        lbb = self.system.lbboundaries

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape2))

        self.check_boundary_flags([1, 0, 2])

    def test_union(self):
        union = espressomd.shapes.Union()
        union.add([self.wall_shape1, self.wall_shape2])
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=union))
        self.check_boundary_flags([1, 0, 1])


@utx.skipIfMissingFeatures(["LB_BOUNDARIES"])
class LBBoundariesCPU(ut.TestCase, LBBoundariesBase):
    lbf = None

    def setUp(self):
        if not self.lbf:
            self.lbf = espressomd.lb.LBFluid(
                visc=1.0,
                dens=1.0,
                agrid=0.5,
                tau=1.0)

        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.lbboundaries.clear()
        self.system.actors.remove(self.lbf)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["LB_BOUNDARIES_GPU"])
class LBBoundariesGPU(ut.TestCase, LBBoundariesBase):
    lbf = None

    def setUp(self):
        if not self.lbf:
            self.lbf = espressomd.lb.LBFluidGPU(
                visc=1.0,
                dens=1.0,
                agrid=0.5,
                tau=1.0)

        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.lbboundaries.clear()
        self.system.actors.remove(self.lbf)


if __name__ == "__main__":
    ut.main()
