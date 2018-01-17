from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
import itertools
import collections


class RandomPairTest(ut.TestCase):
    """This test creates a system of random particles.
       Then the interaction paris for a certain cutoff
       are calculated by brute force in python (pairs_n2),
       and compared to the pairs returned by the cell
       systems, which should be identical. This check is
       repeated for all valid combination of periodicities.

    """
    s = espressomd.System()

    def setUp(self):
        s = self.s
        s.time_step = 1.
        s.box_l = 3 * [10.]
        s.cell_system.skin = 0.0
        s.min_global_cut = 1.5
        n_part = 500

        np.random.seed(2)

        s.part.add(pos=s.box_l * np.random.random((n_part, 3)))

    def tearDown(self):
        self.s.part.clear()

    def pairs_n2(self, dist):
        parts = self.s.part

        pairs = []
        for i in range(len(parts)):
            for j in range(i + 1, len(parts)):
                if self.s.distance(parts[i], parts[j]) < dist:
                    pairs.append((i, j))

        self.assertTrue(len(pairs))
        return set(pairs)

    def check_duplicates(self, l):
        for e in collections.Counter(l).values():
            self.assertEqual(e, 1)

    def check_pairs(self, n2_pairs):
        cs_pairs = self.s.cell_system.get_pairs_(1.5)
        self.check_duplicates(cs_pairs)
        self.assertTrue(len(cs_pairs))
        self.assertEqual(n2_pairs ^ set(cs_pairs), set())

    def check_dd(self, n2_pairs):
        self.s.cell_system.set_domain_decomposition()
        self.check_pairs(n2_pairs)

    def check_layered(self, n2_pairs):
        self.s.cell_system.set_layered()
        self.check_pairs(n2_pairs)

    def check_n_squared(self, n2_pairs):
        self.s.cell_system.set_n_square()
        self.check_pairs(n2_pairs)

    def test(self):
        if espressomd.has_features("PARTIAL_PERIODIC"):
            periods = [0, 1]
        else:
            periods = [1]

        for periodicity in itertools.product(periods, periods, periods):
            self.s.periodicity = periodicity
            n2_pairs = self.pairs_n2(1.5)

            self.check_dd(n2_pairs)
            self.check_layered(n2_pairs)
            self.check_n_squared(n2_pairs)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
