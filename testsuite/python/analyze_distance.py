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
import sys
import unittest as ut
import numpy as np
import espressomd

BOX_L = 50.


@ut.skipIf(not espressomd.has_features("LENNARD_JONES"), "Skipped because LENNARD_JONES turned off.")
class AnalyzeDistance(ut.TestCase):
    system = espressomd.System(box_l=3 * [BOX_L])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(1234)

    def setUp(self):
        self.system.part.add(pos=np.random.random((100, 3)) * BOX_L)

    def tearDown(self):
        self.system.part.clear()

    # python version of the espresso core function
    def min_dist(self):
        r = np.array(self.system.part[:].pos)
        # this generates indices for all i<j combinations
        ij = np.triu_indices(len(r), k=1)
        r_ij = np.fabs(r[ij[0]] - r[ij[1]])
        # check smaller distances via PBC
        r_ij = np.where(
            r_ij > 0.5 * self.system.box_l, self.system.box_l - r_ij, r_ij)
        dist = np.sum(r_ij**2, axis=1)
        return np.sqrt(np.min(dist))

    # python version of the espresso core function
    def nbhood(self, pos, r_catch):
        dist = np.fabs(np.array(self.system.part[:].pos) - pos)
        # check smaller distances via PBC
        dist = np.where(
            dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=1)
        return np.where(dist < r_catch**2)[0]

    # python version of the espresso core function, using pos
    def dist_to_pos(self, pos):
        dist = np.fabs(self.system.part[:].pos - pos)
        # check smaller distances via PBC
        dist = np.where(
            dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=-1)
        return np.sqrt(np.min(dist))

    # python version of the espresso core function, using id
    def dist_to_id(self, id):
        dist = np.fabs(
            np.delete(self.system.part[:].pos, id, axis=0) - self.system.part[id].pos)
        # check smaller distances via PBC
        dist = np.where(
            dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=1)
        return np.sqrt(np.min(dist))

    def test_min_dist(self):
        # try five times
        for i in range(5):
            self.system.part[:].pos = np.random.random(
                (len(self.system.part), 3)) * BOX_L
            self.assertAlmostEqual(self.system.analysis.min_dist(),
                                   self.min_dist(),
                                   delta=1e-7)

    def test_min_dist_empty(self):
        self.system.part.clear()
        self.assertEqual(self.system.analysis.min_dist(), float("inf"))

    def test_nbhood(self):
        # try five times
        for i in range(1, 10, 2):
            self.system.part[:].pos = np.random.random(
                (len(self.system.part), 3)) * BOX_L
            self.assertTrue(
                np.allclose(self.system.analysis.nbhood([i, i, i], i * 2),
                            self.nbhood([i, i, i], i * 2)))

    def test_dist_to_pos(self):
        # try five times
        for i in range(5):
            self.system.part[:].pos = np.random.random(
                (len(self.system.part), 3)) * BOX_L
            self.assertTrue(
                np.allclose(self.system.analysis.dist_to(pos=[i, i, i]),
                            self.dist_to_pos([i, i, i])))

    def test_dist_to_id(self):
        # try five times
        for i in range(5):
            self.system.part[:].pos = np.random.random(
                (len(self.system.part), 3)) * BOX_L
            self.assertAlmostEqual(self.system.analysis.dist_to(id=i),
                                   self.dist_to_id(i))

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
