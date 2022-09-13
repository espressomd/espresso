#
# Copyright (C) 2010-2022 The ESPResSo project
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
import espressomd
import numpy as np
import scipy.spatial
import itertools
import collections
import tests_common


class RandomPairTest(ut.TestCase):

    """This test creates a system of random particles.
       Then the interaction pairs for a certain cutoff
       are calculated by brute force in python (pairs_n2),
       and compared to the pairs returned by the cell
       systems, which should be identical. This check is
       repeated for all valid combinations of periodicities.

    """
    system = espressomd.System(box_l=[10., 15., 15.])

    def setUp(self):
        system = self.system
        system.time_step = .1
        system.cell_system.skin = 0.0
        system.min_global_cut = 1.5
        n_part = 400

        np.random.seed(2)

        positions = system.box_l * np.random.random((n_part, 3))
        system.part.add(pos=positions)
        self.all_pairs = []

        def euclidean_pbc(a, b, box_l):
            vec = np.fmod(a - b + box_l, box_l)
            for i in range(3):
                if vec[i] > box_l[i] / 2.:
                    vec[i] -= box_l[i]
            return np.linalg.norm(vec)

        dist_mat = scipy.spatial.distance.cdist(
            positions, positions, metric=euclidean_pbc,
            box_l=np.copy(system.box_l))
        diag_mask = np.logical_not(np.eye(n_part, dtype=bool))
        dist_crit = (dist_mat < 1.5) * diag_mask
        self.all_pairs = list(zip(*np.nonzero(np.triu(dist_crit))))
        self.assertGreater(len(self.all_pairs), 0)

    def tearDown(self):
        self.system.part.clear()

    def pairs_n2(self, dist):
        # Go through list of all possible pairs for full periodicity
        # and skip those that are not within the desired distance
        # for the current periodicity

        pairs = []
        parts = self.system.part
        for p in self.all_pairs:
            if self.system.distance(parts.by_id(
                    p[0]), parts.by_id(p[1])) <= dist:
                pairs.append(p)
        return set(pairs)

    def check_duplicates(self, l):
        for e in collections.Counter(l).values():
            self.assertEqual(e, 1)

    def check_pairs(self, n2_pairs):
        cs_pairs = self.system.cell_system.get_pairs(1.5)
        self.check_duplicates(cs_pairs)
        self.assertGreater(len(cs_pairs), 0)
        self.assertEqual(n2_pairs ^ set(cs_pairs), set())

    def check_dd(self, n2_pairs):
        self.system.cell_system.set_regular_decomposition()
        self.check_pairs(n2_pairs)

    def check_n_squared(self, n2_pairs):
        self.system.cell_system.set_n_square()
        self.check_pairs(n2_pairs)

    def test(self):
        periods = [0, 1]
        self.system.periodicity = [True, True, True]
        tests_common.check_non_bonded_loop_trace(self, self.system)

        for periodicity in itertools.product(periods, periods, periods):
            self.system.periodicity = periodicity
            n2_pairs = self.pairs_n2(1.5)

            self.check_dd(n2_pairs)
            self.check_n_squared(n2_pairs)


if __name__ == '__main__':
    ut.main()
