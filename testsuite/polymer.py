from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import FeneBond
from espressomd import polymer


class Polymer(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(1234)
    system.set_random_state_PRNG()

    num_mono = 5

    @classmethod
    def setUpClass(self):
        box_l = 20.0
        # start with a small bo
        self.system.box_l = np.array([box_l, box_l, box_l])
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        fene = FeneBond(k=30, d_r_max=2)
        self.fene = fene
        self.system.bonded_inter.add(fene)

    def test(self):
        num_poly = 2
        num_mono = 5
        polymer.create_polymer(start_pos=[1, 1, 1], N_P=num_poly,
                               bond_length=0.9, bond=self.fene,
                               MPC=num_mono,
                               start_id=2)

        # Was the start id considered
        # bond=fene,start_id=2)
        for i in 0, 1:
            self.assertTrue(not self.system.part.exists(i))
        # Were all other particles placed in the correct order
        for i in range(2, 2 + num_mono * num_poly):
            self.assertTrue(self.system.part.exists(i))
        # Total number of particles
        self.assertEqual(len(self.system.part), num_mono * num_poly)

        # Start position
        np.testing.assert_allclose(
            np.copy(self.system.part[2].pos), [1., 1., 1.])

        # Distance between consecutive particles
        for i in range(num_poly):
            first_particle = 2 + num_mono * i
            for j in range(first_particle + 1, first_particle + num_mono):
                print(first_particle, j)
                self.assertAlmostEqual(self.system.distance(
                    self.system.part[j], self.system.part[j - 1]), 0.9, places=5)
        # Test polymer with specified pos2
        self.system.part.clear()
        polymer.create_polymer(start_pos=[1, 1, 1], pos2=[1.9, 1, 1], N_P=1,
                               bond_length=0.9, bond=self.fene,
                               MPC=num_mono,
                               start_id=2, angle2=1)

        np.testing.assert_allclose(
            np.copy(self.system.part[2].pos), [1., 1., 1.])
        np.testing.assert_allclose(
            np.copy(self.system.part[3].pos), [1.9, 1., 1.])


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
