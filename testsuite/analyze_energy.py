from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import HarmonicBond

@ut.skipIf(not espressomd.has_features("LENNARD_JONES"), "Skipped because LENNARD_JONES turned off.")
class AnalyzeEnergy(ut.TestCase):
    system = espressomd.System()
    harmonic = HarmonicBond(r_0=0.0, k=3)

    @classmethod
    def setUpClass(self):
        box_l = 10.0
        self.system.box_l = [box_l, box_l, box_l]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.non_bonded_inter[1, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.part.add(id=0, pos=[1, 2, 2], type=0)
        self.system.part.add(id=1, pos=[5, 2, 2], type=0)
        self.system.thermostat.set_langevin(kT=0., gamma=1.)
        self.system.bonded_inter.add(self.harmonic)

    def test_kinetic(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [5, 2, 2]
        self.system.part[0].v = [3, 4, 5]
        self.system.part[1].v = [0, 0, 0]
        # initial, single moving particle
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 25.))
        self.assertTrue(np.allclose(energy["kinetic"], 25.))
        self.assertTrue(np.allclose(energy["bonded"], 0.))
        self.assertTrue(np.allclose(energy["non_bonded"], 0.))
        # initial, two moving particles
        self.system.part[1].v = [3, 4, 5]
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 50.))
        self.assertTrue(np.allclose(energy["kinetic"], 50.))
        self.assertTrue(np.allclose(energy["bonded"], 0.))
        self.assertTrue(np.allclose(energy["non_bonded"], 0.))
        # rest
        self.system.integrator.run(10000)
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 0.))
        self.assertTrue(np.allclose(energy["kinetic"], 0.))
        self.assertTrue(np.allclose(energy["bonded"], 0.))
        self.assertTrue(np.allclose(energy["non_bonded"], 0.))

    def test_non_bonded(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [2, 2, 2]
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 1.))
        self.assertTrue(np.allclose(energy["kinetic"], 0.))
        self.assertTrue(np.allclose(energy["bonded"], 0.))
        self.assertTrue(np.allclose(energy["non_bonded"], 1.))
        # add another pair of particles
        self.system.part.add(id=3, pos=[1, 5, 5], type=1)
        self.system.part.add(id=4, pos=[2, 5, 5], type=1)
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 2.))
        self.assertTrue(np.allclose(energy["kinetic"], 0.))
        self.assertTrue(np.allclose(energy["bonded"], 0.))
        self.assertTrue(np.allclose(energy["non_bonded"], 2.))
        self.system.part[3].remove()
        self.system.part[4].remove()

    def test_bonded(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [3, 2, 2]
        self.system.part[0].v = [0, 0, 0]
        self.system.part[1].v = [0, 0, 0]
        # single bond
        self.system.part[0].add_bond((self.harmonic, 1))
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 6))
        self.assertTrue(np.allclose(energy["kinetic"], 0.))
        self.assertTrue(np.allclose(energy["bonded"], 6))
        self.assertTrue(np.allclose(energy["non_bonded"], 0.))
        # two bonds
        self.system.part[0].add_bond((self.harmonic, 1))
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 12))
        self.assertTrue(np.allclose(energy["kinetic"], 0.))
        self.assertTrue(np.allclose(energy["bonded"], 12))
        self.assertTrue(np.allclose(energy["non_bonded"], 0.))
        # bonds deleted
        self.system.part[0].delete_all_bonds()
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 0.))
        self.assertTrue(np.allclose(energy["kinetic"], 0.))
        self.assertTrue(np.allclose(energy["bonded"], 0))
        self.assertTrue(np.allclose(energy["non_bonded"], 0.))

    def test_all(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [2, 2, 2]
        self.system.part[0].v = [3, 4, 5]
        self.system.part[1].v = [3, 4, 5]
        # single bond
        self.system.part[0].add_bond((self.harmonic, 1))
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 50. + 3. / 2. + 1.))
        self.assertTrue(np.allclose(energy["kinetic"], 50.))
        self.assertTrue(np.allclose(energy["bonded"], 3. / 2.))
        self.assertTrue(np.allclose(energy["non_bonded"], 1.))
        # two bonds
        self.system.part[0].add_bond((self.harmonic, 1))
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 50. + 3 + 1.))
        self.assertTrue(np.allclose(energy["kinetic"], 50.))
        self.assertTrue(np.allclose(energy["bonded"], 3.))
        self.assertTrue(np.allclose(energy["non_bonded"], 1.))
        # add another pair of particles
        self.system.part.add(id=3, pos=[1, 5, 5], type=1)
        self.system.part.add(id=4, pos=[2, 5, 5], type=1)
        energy = self.system.analysis.energy()
        self.assertTrue(np.allclose(energy["total"], 50. + 3 + (1. + 1.)))
        self.assertTrue(np.allclose(energy["kinetic"], 50.))
        self.assertTrue(np.allclose(energy["bonded"], 3.))
        self.assertTrue(np.allclose(energy["non_bonded"], 1. + 1.))
        self.system.part[3].remove()
        self.system.part[4].remove()
        self.system.part[0].delete_all_bonds()

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
