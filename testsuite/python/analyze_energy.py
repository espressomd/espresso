# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.interactions
import espressomd.electrostatics


@utx.skipIfMissingFeatures("LENNARD_JONES")
class AnalyzeEnergy(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    harmonic = espressomd.interactions.HarmonicBond(r_0=0.0, k=3)

    @classmethod
    def setUpClass(cls):
        box_l = 20
        cls.system.box_l = [box_l, box_l, box_l]
        cls.system.cell_system.skin = 0.4
        cls.system.time_step = 0.01
        cls.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        cls.system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        cls.system.non_bonded_inter[1, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        cls.system.thermostat.set_langevin(kT=0., gamma=1., seed=42)
        cls.system.bonded_inter[5] = cls.harmonic

    def setUp(self):
        self.system.part.add(pos=[1, 2, 2], type=0)
        self.system.part.add(pos=[5, 2, 2], type=0)

    def tearDown(self):
        self.system.part.clear()

    def test_kinetic(self):
        p0, p1 = self.system.part[:]
        p0.pos = [1, 2, 2]
        p1.pos = [5, 2, 2]
        # single moving particle
        p0.v = [3, 4, 5]
        p1.v = [0, 0, 0]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # two moving particles
        p1.v = [3, 4, 5]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)

    def test_non_bonded(self):
        p0, p1 = self.system.part[:]
        p0.pos = [1, 2, 2]
        p1.pos = [2, 2, 2]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 1., delta=1e-5)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        # add another pair of particles
        self.system.part.add(pos=[3, 2, 2], type=1)
        self.system.part.add(pos=[4, 2, 2], type=1)
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded", 0, 1], 1., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded", 0, 0]
                               + energy["non_bonded", 0, 1]
                               + energy["non_bonded", 1, 1], energy["total"], delta=1e-7)

    def test_bonded(self):
        p0, p1 = self.system.part[:]
        p0.pos = [1, 2, 2]
        p1.pos = [3, 2, 2]
        # single bond
        p0.add_bond((self.harmonic, p1))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 6, delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 6, delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # two bonds
        p1.add_bond((self.harmonic, p0))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 12, delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 12, delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # bonds deleted
        p0.delete_all_bonds()
        p1.delete_all_bonds()
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0, delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)

    def test_all(self):
        p0, p1 = self.system.part[:]
        p0.pos = [1, 2, 2]
        p1.pos = [2, 2, 2]
        p0.v = [3, 4, 5]
        p1.v = [3, 4, 5]
        # single bond
        p0.add_bond((self.harmonic, p1))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50. + 3. / 2. + 1., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3. / 2., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        # two bonds
        p1.add_bond((self.harmonic, p0))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50. + 3 + 1., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        # add another pair of particles
        self.system.part.add(pos=[1, 5, 5], type=1)
        self.system.part.add(pos=[2, 5, 5], type=1)
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(
            energy["total"], 50. + 3 + (1. + 1.), delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1. + 1., delta=1e-7)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS", "P3M"])
    def test_electrostatics(self):
        p0, p1 = self.system.part[:]
        p0.pos = [1, 2, 2]
        p1.pos = [3, 2, 2]
        p0.q = 1
        p1.q = -1
        p3m = espressomd.electrostatics.P3M(
            prefactor=1.0,
            accuracy=9.910945054074526e-08,
            mesh=[22, 22, 22],
            cao=7,
            r_cut=8.906249999999998,
            alpha=0.387611049779351,
            tune=False)
        self.system.actors.add(p3m)

        # did not verify if this is correct, but looks pretty good (close to
        # 1/2)
        u_p3m = -0.501062398379
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], u_p3m, delta=1e-5)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0, delta=1e-7)
        self.assertAlmostEqual(energy["coulomb"], u_p3m, delta=1e-5)


if __name__ == "__main__":
    ut.main()
