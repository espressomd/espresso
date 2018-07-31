from __future__ import print_function
import sys
import unittest as ut
import espressomd
from espressomd.interactions import HarmonicBond


@ut.skipIf(not espressomd.has_features("LENNARD_JONES"), "Skipped because LENNARD_JONES turned off.")
class AnalyzeEnergy(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    harmonic = HarmonicBond(r_0=0.0, k=3)

    @classmethod
    def setUpClass(self):
        box_l = 20
        self.system.box_l = [box_l, box_l, box_l]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.non_bonded_inter[1, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.thermostat.set_langevin(kT=0., gamma=1.)
        self.system.bonded_inter.add(self.harmonic)

    def setUp(self):
        self.system.part.clear()
        self.system.part.add(id=0, pos=[1, 2, 2], type=0)
        self.system.part.add(id=1, pos=[5, 2, 2], type=0)

    def test_kinetic(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [5, 2, 2]
        self.system.part[0].v = [3, 4, 5]
        self.system.part[1].v = [0, 0, 0]
        # single moving particle
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # two moving particles
        self.system.part[1].v = [3, 4, 5]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        self.system.part[0].v = [0, 0, 0]
        self.system.part[1].v = [0, 0, 0]


    def test_non_bonded(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [2, 2, 2]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 1., delta=1e-5)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        # add another pair of particles
        self.system.part.add(id=2, pos=[3, 2, 2], type=1)
        self.system.part.add(id=3, pos=[4, 2, 2], type=1)
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded", 0, 1], energy["non_bonded", 1, 0], delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded", 0, 0]
                               + energy["non_bonded", 0, 1]
                               + energy["non_bonded", 1, 1], energy["total"], delta=1e-7)
        self.system.part[2].remove()
        self.system.part[3].remove()

    def test_bonded(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [3, 2, 2]
        self.system.part[0].v = [0, 0, 0]
        self.system.part[1].v = [0, 0, 0]
        # single bond
        self.system.part[0].add_bond((self.harmonic, 1))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 6, delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 6, delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # two bonds
        self.system.part[1].add_bond((self.harmonic, 0))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 12, delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 12, delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # bonds deleted
        self.system.part[0].delete_all_bonds()
        self.system.part[1].delete_all_bonds()
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0, delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)

    def test_all(self):
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [2, 2, 2]
        self.system.part[0].v = [3, 4, 5]
        self.system.part[1].v = [3, 4, 5]
        # single bond
        self.system.part[0].add_bond((self.harmonic, 1))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50. + 3. / 2. + 1., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3. / 2., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        # two bonds
        self.system.part[1].add_bond((self.harmonic, 0))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50. + 3 + 1., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        # add another pair of particles
        self.system.part.add(id=2, pos=[1, 5, 5], type=1)
        self.system.part.add(id=3, pos=[2, 5, 5], type=1)
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(
            energy["total"], 50. + 3 + (1. + 1.), delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1. + 1., delta=1e-7)
        self.system.part[2].remove()
        self.system.part[3].remove()
        self.system.part[0].delete_all_bonds()

    features = ["ELECTROSTATICS", "P3M"]

    @ut.skipIf(not espressomd.has_features(*features),
               "Missing features: " + ", ".join(espressomd.missing_features(*features)))
    def test_electrostatics(self):

        from espressomd import electrostatics

        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [3, 2, 2]
        self.system.part[0].q = 1
        self.system.part[1].q = -1
        p3m = electrostatics.P3M(prefactor=1.0,
                                 accuracy=9.910945054074526e-08,
                                 mesh=[22, 22, 22],
                                 cao=7,
                                 r_cut=8.906249999999998,
                                 alpha=0.387611049779351,
                                 tune=False)
        self.system.actors.add(p3m)

        # did not verify if this is correct, but looks pretty good (close to 1/2)
        u_p3m = -0.501062398379
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], u_p3m, delta=1e-5)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0, delta=1e-7)
        self.assertAlmostEqual(energy["coulomb"], u_p3m, delta=1e-7)
        self.system.part[0].q = 0
        self.system.part[1].q = 0
        self.system.part[0].pos = [1, 2, 2]
        self.system.part[1].pos = [5, 2, 2]


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
