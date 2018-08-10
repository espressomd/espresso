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
    system = espressomd.System(box_l = 3 * [10.])
        
    def setUp(self):
        s = self.system
        s.time_step = .1
        s.cell_system.skin = 0.0
        s.min_global_cut = 1.5
        n_part = 500

        np.random.seed(2)

        s.part.add(pos=s.box_l * np.random.random((n_part, 3)))
        self.all_pairs=[]

        dist_func=self.system.distance
        for pair in self.system.part.pairs():
                if dist_func(pair[0], pair[1]) < 1.5:
                    self.all_pairs.append((pair[0].id,pair[1].id))

        self.all_pairs=set(self.all_pairs)
        self.assertTrue(len(self.all_pairs))



    def tearDown(self):
        self.system.part.clear()

    def pairs_n2(self, dist):
        # Go through list of all possible pairs for full periodicy
        # and skip those that ar not within the desired distance
        # for the current periodicity

        pairs=[]
        parts=self.system.part
        for p in self.all_pairs:
            if self.system.distance(parts[p[0]],parts[p[1]]) <=dist:
                pairs.append(p)
        return set(pairs)

    def check_duplicates(self, l):
        for e in collections.Counter(l).values():
            self.assertEqual(e, 1)

    def check_pairs(self, n2_pairs):
        cs_pairs = self.system.cell_system.get_pairs_(1.5)
        self.check_duplicates(cs_pairs)
        self.assertTrue(len(cs_pairs))
        self.assertEqual(n2_pairs ^ set(cs_pairs), set())

    def check_dd(self, n2_pairs):
        self.system.cell_system.set_domain_decomposition()
        self.check_pairs(n2_pairs)

    def check_layered(self, n2_pairs):
        self.system.cell_system.set_layered()
        self.check_pairs(n2_pairs)

    def check_n_squared(self, n2_pairs):
        self.system.cell_system.set_n_square()
        self.check_pairs(n2_pairs)

    def test(self):
        if espressomd.has_features("PARTIAL_PERIODIC"):
            periods = [0, 1]
        else:
            periods = [1]

        for periodicity in itertools.product(periods, periods, periods):
            self.system.periodicity = periodicity
            n2_pairs = self.pairs_n2(1.5)

            self.check_dd(n2_pairs)
            self.check_layered(n2_pairs)
            self.check_n_squared(n2_pairs)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
