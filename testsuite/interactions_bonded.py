#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from __future__ import print_function

import unittest as ut
import numpy as np

import espressomd
import tests_common


class InteractionsBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(seed=system.seed)

    box_l = 10.

    start_pos = np.random.rand(3) * box_l
    axis = np.random.rand(3)
    axis /= np.linalg.norm(axis)
    step = axis * 0.01
    step_width = np.linalg.norm(step)

    def setUp(self):

        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos, type=0)

    def tearDown(self):

        self.system.part.clear()

    # Test Harmonic Bond
    def test_harmonic(self):

        hb_k = 5
        hb_r_0 = 1.5
        hb_r_cut = 3.355

        hb = espressomd.interactions.HarmonicBond(
            k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)
        self.system.bonded_inter.add(hb)
        self.system.part[0].add_bond((hb, 1))

        for i in range(335):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = tests_common.harmonic_potential(
                scalar_r=(i + 1) * self.step_width, k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.harmonic_force(scalar_r=(i + 1) * self.step_width,
                               k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)

            # Check that energies match, ...
            np.testing.assert_almost_equal(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            f1_sim_copy = np.copy(f1_sim)
            np.testing.assert_almost_equal(f1_sim_copy, f1_ref)

        # Check that bond breaks when distance > r_cut
        self.system.part[1].pos = self.system.part[1].pos + self.step
        with self.assertRaisesRegexp(Exception, "Encoutered errors during integrate"):
            self.system.integrator.run(recalc_forces=True, steps=0)

    # Test Fene Bond
    def test_fene(self):
        fene_k = 23.15
        fene_d_r_max = 3.355
        fene_r_0 = 1.1

        fene = espressomd.interactions.FeneBond(
            k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)
        self.system.bonded_inter.add(fene)
        self.system.part[0].add_bond((fene, 1))

        for i in range(445):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = tests_common.fene_potential(
                scalar_r=(i + 1) * self.step_width, k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.fene_force(scalar_r=(i + 1) * self.step_width,
                           k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)

            # Check that energies match, ...
            np.testing.assert_almost_equal(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            f1_sim_copy = np.copy(f1_sim)
            np.testing.assert_almost_equal(f1_sim_copy, f1_ref, decimal=5)

        # Check that bond breaks when distance > r_cut
        self.system.part[1].pos = self.system.part[1].pos + self.step
        with self.assertRaisesRegexp(Exception, "Encoutered errors during integrate"):
            self.system.integrator.run(recalc_forces=True, steps=0)

    @ut.skipIf(not espressomd.has_features(["ELECTROSTATICS"]),
               "ELECTROSTATICS feature is not available, skipping coulomb test.")
    def test_coulomb(self):
        coulomb_k = 1
        q1 = 1
        q2 = -1
        self.system.part[0].q = q1
        self.system.part[1].q = q2

        coulomb = espressomd.interactions.BondedCoulomb(prefactor=coulomb_k)
        self.system.bonded_inter.add(coulomb)
        self.system.part[0].add_bond((coulomb, 1))

        for i in range(445):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = tests_common.coulomb_potential(
                scalar_r=(i + 1) * self.step_width, k=coulomb_k, q1=q1, q2=q2)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * tests_common.coulomb_force(
                scalar_r=(i + 1) * self.step_width,
                                k=coulomb_k, q1=q1, q2=q2)

            # Check that energies match, ...
            np.testing.assert_almost_equal(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            f1_sim_copy = np.copy(f1_sim)
            np.testing.assert_almost_equal(f1_sim_copy, f1_ref)

if __name__ == '__main__':
    ut.main()
