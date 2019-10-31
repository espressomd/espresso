#
# Copyright (C) 2013-2019 The ESPResSo project
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
import unittest_decorators as utx
import espressomd
from espressomd.interactions import FeneBond
from espressomd import pair_criteria


class PairCriteria(ut.TestCase):

    """Tests interface and implementation of pair criteria"""

    es = espressomd.System(box_l=[1., 1., 1.])

    f1 = FeneBond(k=1, d_r_max=0.05)
    es.bonded_inter.add(f1)
    f2 = FeneBond(k=1, d_r_max=0.05)
    es.bonded_inter.add(f2)
    es.part.add(id=0, pos=(0, 0, 0))
    es.part.add(id=1, pos=(0.91, 0, 0))
    p1 = es.part[0]
    p2 = es.part[1]
    epsilon = 1E-8

    def test_distance_crit_periodic(self):
        dc = pair_criteria.DistanceCriterion(cut_off=0.1)
        # Interface
        self.assertEqual(list(dc.get_params().keys()), ["cut_off", ])
        self.assertTrue(abs(dc.get_params()["cut_off"] - 0.1) < self.epsilon)

        # Decisions
        # Periodic system. Particles in range via minimum image convention
        self.es.periodicity = (1, 1, 1)
        self.assertTrue(dc.decide(self.p1, self.p2))
        self.assertTrue(dc.decide(self.p1.id, self.p2.id))

    def test_distance_crit_non_periodic(self):
        dc = pair_criteria.DistanceCriterion(cut_off=0.1)

        # Non-periodic system. Particles out of range
        self.es.periodicity = (0, 0, 0)
        self.assertTrue(not dc.decide(self.p1, self.p2))
        self.assertTrue(not dc.decide(self.p1.id, self.p2.id))

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_energy_crit(self):
        # Setup purely repulsive lj
        self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=0.11, epsilon=1, cutoff=2**(1. / 6.) * 0.11, shift="auto")
        ec = pair_criteria.EnergyCriterion(cut_off=0.001)
        # Interface
        self.assertEqual(list(ec.get_params().keys()), ["cut_off", ])
        self.assertTrue(abs(ec.get_params()["cut_off"] - 0.001) < self.epsilon)

        # Decisions
        # Periodic system. Particles in range via minimum image convention
        self.es.periodicity = (1, 1, 1)
        self.assertTrue(ec.decide(self.p1, self.p2))
        self.assertTrue(ec.decide(self.p1.id, self.p2.id))

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_energy_crit_non_periodic(self):
        # Setup purely repulsive lj
        self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=0.11, epsilon=1, cutoff=2**(1. / 6.) * 0.11, shift="auto")
        ec = pair_criteria.EnergyCriterion(cut_off=0.001)
        # Interface
        self.assertEqual(list(ec.get_params().keys()), ["cut_off", ])
        self.assertTrue(abs(ec.get_params()["cut_off"] - 0.001) < self.epsilon)

        # Non-periodic system. Particles out of range
        self.es.periodicity = (0, 0, 0)
        self.assertTrue(not ec.decide(self.p1, self.p2))
        self.assertTrue(not ec.decide(self.p1.id, self.p2.id))

    def test_bond_crit(self):
        bc = pair_criteria.BondCriterion(bond_type=0)
        # Interface
        self.assertEqual(list(bc.get_params().keys()), ["bond_type", ])
        self.assertEqual(bc.get_params()["bond_type"], 0)

        # Decisions
        # No bond yet. Should return false
        self.assertTrue(not bc.decide(self.p1, self.p2))
        self.assertTrue(not bc.decide(self.p1.id, self.p2.id))

        # Add bond. Then the criterion should match
        self.es.part[0].bonds = ((0, 1),)
        self.assertTrue(bc.decide(self.p1, self.p2))
        self.assertTrue(bc.decide(self.p1.id, self.p2.id))

        # Place bond on the 2nd particle. The criterion should still match
        self.es.part[0].bonds = ()
        self.es.part[1].bonds = ((0, 0),)
        self.assertTrue(bc.decide(self.p1, self.p2))
        self.assertTrue(bc.decide(self.p1.id, self.p2.id))


if __name__ == "__main__":
    ut.main()
