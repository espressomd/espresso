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

    @classmethod
    def setUpClass(self):
        s = self.system
        s.time_step = 1.
        s.cell_system.skin = 0.0
        s.box_l = 10. * self.system.cell_system.get_state()['node_grid']
        s.min_global_cut = 1.5
        n_part = 100 * s.cell_system.get_state()['n_nodes']

        print("n_part", n_part)

        np.random.seed(2)

        s.part.add(pos=s.box_l * np.random.random((n_part, 3)))

    def reset_parts(self):
        s = self.system

        n_part = len(s.part)
        s.part[:].pos = s.box_l * np.random.random((n_part, 3))
        s.part[:].v = (-0.5 + np.random.random((n_part, 3)))

    def pairs_n2(self, dist):
        parts = self.system.part

        pairs = []
        for i in range(len(parts)):
            for j in range(i + 1, len(parts)):
                if self.system.distance(parts[i], parts[j]) < dist:
                    pairs.append((i, j))

        return set(pairs)

    def check_duplicates(self, l):
        for e in collections.Counter(l).values():
            self.assertEqual(e, 1)

    def check_pairs(self):
        n2_pairs = self.pairs_n2(1.5)
        cs_pairs = self.system.cell_system.get_pairs_(1.5)
        self.check_duplicates(cs_pairs)
        self.assertEqual(n2_pairs, set(cs_pairs))

    def run_checks(self):
        if espressomd.has_features("PARTIAL_PERIODIC"):
            periods = [0, 1]
        else:
            periods = [1]

        for periodicity in itertools.product(periods, periods, periods):
            print(periodicity)
            self.system.periodicity = periodicity
            self.reset_parts()

            self.check_pairs()

            self.system.integrator.run(1000)

            self.check_pairs()

    def test_dd(self):
        self.system.box_l = 10. * self.system.cell_system.get_state()['node_grid']
        self.system.cell_system.set_domain_decomposition()
        print(self.system.cell_system.get_state())
        self.run_checks()

    def test_layered(self):
        self.system.box_l = [10., 10, 10. * self.system.cell_system.get_state()['n_nodes']]
        self.system.cell_system.set_layered()
        print(self.system.cell_system.get_state())
        self.run_checks()

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
