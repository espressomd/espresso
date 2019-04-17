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
from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
import espressomd.lb
from espressomd.shapes import Wall
import espressomd.lbboundaries
from itertools import product


class LBBoundariesBase(object):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    wall_shape1 = Wall(normal=[1., 0., 0.], dist=2.5)
    wall_shape2 = Wall(normal=[-1., 0., 0.], dist=-7.5)

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

        self.assertFalse(b1 in lbb)
        self.assertTrue(b2 in lbb)

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

        b1 = lbb.add(
            espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        b2 = lbb.add(
            espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))

        lbb.clear()

        self.assertTrue(lbb.empty())

    def test_boundary_flags(self):
        lbb = self.system.lbboundaries

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape2))

        rng = range(20)
        lbf = self.lbf

        for i in product(range(0, 5), rng, rng):
            self.assertEqual(lbf[i].boundary, 1)

        for i in product(range(5, 15), rng, rng):
            self.assertEqual(lbf[i].boundary, 0)

        for i in product(range(15, 20), rng, rng):
            self.assertEqual(lbf[i].boundary, 2)

        lbb.clear()
        for i in product(rng, rng, rng):
            self.assertEqual(lbf[i].boundary, 0)


@ut.skipIf(not espressomd.has_features(["LB_BOUNDARIES"]),
           "Features not available, skipping test!")
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


@ut.skipIf(
    not espressomd.gpu_available() or not espressomd.has_features(
        ["LB_BOUNDARIES_GPU"]),
           "Features or not available, skipping test!")
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
