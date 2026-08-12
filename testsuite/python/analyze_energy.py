#
# Copyright (C) 2010-2026 The ESPResSo project
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

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.interactions
import espressomd.electrostatics


@utx.skipIfMissingFeatures("LENNARD_JONES")
class AnalyzeEnergy(ut.TestCase):
    system = espressomd.System(box_l=[20.0, 20.0, 20.0])
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    system.thermostat.set_langevin(kT=0., gamma=1., seed=42)
    harmonic = espressomd.interactions.HarmonicBond(r_0=0.0, k=3)
    system.bonded_inter[5] = harmonic

    def setUp(self):
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.non_bonded_inter[1, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.part.add(pos=[1, 2, 2], type=0, mol_id=6)
        self.system.part.add(pos=[5, 2, 2], type=0, mol_id=6)

    def tearDown(self):
        self.system.non_bonded_inter.reset()
        self.system.part.clear()
        if espressomd.has_features(["ELECTROSTATICS"]):
            self.system.electrostatics.clear()

    def test_kinetic(self):
        p0, p1 = self.system.part.all()
        p0.pos = [1, 2, 2]
        p1.pos = [5, 2, 2]
        # single moving particle
        p0.v = [3, 4, 5]
        p1.v = [0, 0, 0]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic_lin"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic_rot"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 25., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        # two moving particles
        p1.v = [3, 4, 5]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic_lin"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic_rot"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0., delta=1e-7)
        if espressomd.has_features(["ROTATION"]):
            p0.omega_lab = [1, 2, 3]
            p1.omega_lab = [1, 2, 3]
            p0.rotation = [True, True, True]
            p1.rotation = [False, False, False]
            energy = self.system.analysis.energy()
            self.assertAlmostEqual(energy["kinetic_lin"], 50., delta=1e-7)
            self.assertAlmostEqual(energy["kinetic_rot"], 7., delta=1e-7)
            self.assertAlmostEqual(energy["kinetic"], 57., delta=1e-7)

    def test_non_bonded(self):
        p0, p1 = self.system.part.all()
        p0.pos = [1, 2, 2]
        p1.pos = [2, 2, 2]
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 1., delta=1e-5)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded_intra"], 1., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded_inter"], 0., delta=1e-7)
        # Test the single particle energy function
        self.assertAlmostEqual(energy["non_bonded"], 0.5 * sum(
            [self.system.analysis.particle_non_bonded_energy(p) for p in self.system.part.all()]), delta=1e-7)
        # add another pair of particles
        self.system.part.add(pos=[3, 2, 2], type=1, mol_id=7)
        self.system.part.add(pos=[4, 2, 2], type=1, mol_id=7)
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded", 0, 1], 1., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded", 0, 0]
                               + energy["non_bonded", 0, 1]
                               + energy["non_bonded", 1, 1], energy["total"], delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded_intra", 0, 0], 1., delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded_intra", 1, 1], 1., delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded_intra", 0, 1], 0., delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded_inter", 0, 0], 0., delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded_inter", 1, 1], 0., delta=1e-7)
        self.assertAlmostEqual(
            energy["non_bonded_inter", 0, 1], 1., delta=1e-7)
        # Test the single particle energy function
        self.assertAlmostEqual(energy["non_bonded"], 0.5 * sum(
            [self.system.analysis.particle_non_bonded_energy(p) for p in self.system.part.all()]), delta=1e-7)

    def test_bonded(self):
        p0, p1 = self.system.part.all()
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
        p0, p1 = self.system.part.all()
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
        self.assertAlmostEqual(energy["non_bonded"], 0.5 * sum(
            [self.system.analysis.particle_non_bonded_energy(p) for p in self.system.part.all()]), delta=1e-7)
        if espressomd.has_features(["VIRTUAL_SITES"]):
            self.assertAlmostEqual(energy["virtual_sites"], 0., delta=1e-7)
        if espressomd.has_features(["DPD"]):
            self.assertAlmostEqual(energy["dpd"], 0., delta=1e-7)
        # two bonds
        p1.add_bond((self.harmonic, p0))
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], 50. + 3 + 1., delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0.5 * sum(
            [self.system.analysis.particle_non_bonded_energy(p) for p in self.system.part.all()]), delta=1e-7)
        # add another pair of particles
        self.system.part.add(pos=[1, 5, 5], type=1)
        self.system.part.add(pos=[2, 5, 5], type=1)
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(
            energy["total"], 50. + 3 + (1. + 1.), delta=1e-7)
        self.assertAlmostEqual(energy["kinetic"], 50., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 3., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 1. + 1., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0.5 * sum(
            [self.system.analysis.particle_non_bonded_energy(p) for p in self.system.part.all()]), delta=1e-7)
        # check effect of particle resort
        p0_energy_old = self.system.analysis.particle_non_bonded_energy(p0)
        p0.pos = p0.pos  # trigger particle resort
        p0_energy_new = self.system.analysis.particle_non_bonded_energy(p0)
        self.assertAlmostEqual(p0_energy_new, p0_energy_old, delta=1e-7)

    def check_electrostatics(self, gpu):
        p0, p1 = self.system.part.all()
        p0.pos = [1, 2, 2]
        p1.pos = [3, 2, 2]
        p0.q = 1
        p1.q = -1
        p3m = espressomd.electrostatics.P3M(
            prefactor=1.0,
            accuracy=1e-7,
            mesh=[22, 22, 22],
            cao=7,
            gpu=gpu,
            r_cut=8.90625,
            alpha=0.38761105,
            tune=False)
        self.system.electrostatics.solver = p3m

        # did not verify if this is correct, but looks pretty good (close to
        # 1/2)
        u_p3m = -0.501062398379
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["total"], u_p3m, delta=1e-5)
        self.assertAlmostEqual(energy["kinetic"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["bonded"], 0., delta=1e-7)
        self.assertAlmostEqual(energy["non_bonded"], 0, delta=1e-7)
        self.assertAlmostEqual(energy["coulomb"], u_p3m, delta=1e-5)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_electrostatics_cpu(self):
        self.check_electrostatics(False)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    def test_electrostatics_gpu(self):
        self.check_electrostatics(True)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_particle_energy(self):
        self.system.non_bonded_inter.reset()
        self.system.part.clear()
        get_non_bonded_energy = self.system.analysis.particle_non_bonded_energy
        tol = 1e-5
        p1 = self.system.part.add(pos=[0., 0., 0.], type=0, q=+1.)
        p2 = self.system.part.add(pos=[1., 0., 0.], type=0, q=-1.)
        self.assertEqual(get_non_bonded_energy(p1), 0.)
        self.assertEqual(get_non_bonded_energy(p2), 0.)
        # check short-range electrostatics energy is excluded
        self.system.electrostatics.solver = espressomd.electrostatics.DH(
            prefactor=1., kappa=1., r_cut=2.)
        coulomb_energy = self.system.analysis.energy()["coulomb"]
        self.assertAlmostEqual(coulomb_energy, -np.exp(-1.), delta=tol)
        self.assertEqual(get_non_bonded_energy(p1), 0.)
        self.assertEqual(get_non_bonded_energy(p2), 0.)
        if espressomd.has_features(["THOLE"]):
            # check Thole correction is excluded, despite being a non-bonded IA
            self.system.non_bonded_inter[0, 0].thole.set_params(
                scaling_coeff=2., q1q2=p1.q * p2.q)
            coulomb_energy = self.system.analysis.energy()["coulomb"]
            nbonded_energy = self.system.analysis.energy()["non_bonded"]
            self.assertAlmostEqual(coulomb_energy, -np.exp(-1.), delta=tol)
            self.assertAlmostEqual(nbonded_energy, 2. * np.exp(-3.), delta=tol)
            self.assertEqual(get_non_bonded_energy(p1), 0.)
            self.assertEqual(get_non_bonded_energy(p2), 0.)

    def test_fallthrough(self):
        self.assertIsNone(
            self.system.analysis.call_method("_observable_stat_test_fallthrough"))


if __name__ == "__main__":
    ut.main()
