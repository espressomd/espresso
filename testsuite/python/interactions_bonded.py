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
import numpy as np

import espressomd
import espressomd.electrostatics
import tests_common


class InteractionsBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[17.0, 9.0, 8.0])
    np.random.seed(seed=42)

    box_l = 10.

    start_pos = np.random.rand(3) * box_l
    axis = np.random.rand(3)
    axis /= np.linalg.norm(axis)
    steps = 10

    def setUp(self):
        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = .2

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
        self.run_test(hb,
                      lambda r: tests_common.harmonic_force(
                          scalar_r=r, k=hb_k, r_0=hb_r_0),
                      lambda r: tests_common.harmonic_potential(
                          scalar_r=r, k=hb_k, r_0=hb_r_0),
                      0.01, hb_r_cut, True)

    # Test Fene Bond
    def test_fene(self):
        fene_k = 23.15
        fene_d_r_max = 3.355
        fene_r_0 = 1.1

        fene = espressomd.interactions.FeneBond(
            k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)
        self.run_test(fene,
                      lambda r: tests_common.fene_force(
                          scalar_r=r, k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0),
                      lambda r: tests_common.fene_potential(
                          scalar_r=r, k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0),
                      0.01, fene_r_0 + fene_d_r_max, True)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_coulomb(self):
        coulomb_k = 1
        q1 = 1
        q2 = -1
        self.system.part[0].q = q1
        self.system.part[1].q = q2
        self.run_test(
            espressomd.interactions.BondedCoulomb(prefactor=coulomb_k),
            lambda r: tests_common.coulomb_force(r, coulomb_k, q1, q2),
            lambda r: tests_common.coulomb_potential(r, coulomb_k, q1, q2),
            0.01, self.system.box_l[0] / 3)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_coulomb_sr(self):
        # with negated actual charges and only short range int: cancels out all
        # interactions
        q1 = 1.2
        q2 = -q1
        self.system.part[0].q = q1
        self.system.part[1].q = q2
        r_cut = 2

        sr_solver = espressomd.electrostatics.DH(
            prefactor=2, kappa=0.8, r_cut=r_cut)
        self.system.actors.add(sr_solver)
        coulomb_sr = espressomd.interactions.BondedCoulombSRBond(
            q1q2=- q1 * q2)

        # no break test, bond can't break. it extends as far as the short range
        # part of the electrostatics actor
        self.run_test(
            coulomb_sr,
            lambda r: [0., 0., 0.],
            lambda r: 0,
            0.01,
            r_cut,
            test_breakage=False)

    def test_quartic(self):
        """Tests the Quartic bonded interaction by comparing the potential and
           force against the analytic values"""

        quartic_k0 = 2.
        quartic_k1 = 5.
        quartic_r = 0.5
        quartic_r_cut = self.system.box_l[0] / 3.

        quartic = espressomd.interactions.QuarticBond(k0=quartic_k0,
                                                      k1=quartic_k1,
                                                      r=quartic_r,
                                                      r_cut=quartic_r_cut)

        self.run_test(quartic,
                      lambda r: tests_common.quartic_force(
                          k0=quartic_k0, k1=quartic_k1, r=quartic_r, r_cut=quartic_r_cut, scalar_r=r),
                      lambda r: tests_common.quartic_potential(
                          k0=quartic_k0, k1=quartic_k1, r=quartic_r, r_cut=quartic_r_cut, scalar_r=r),
                      0.01, quartic_r_cut, True)

    def run_test(self, bond_instance, force_func, energy_func, min_dist,
                 cutoff, test_breakage=False):
        self.system.bonded_inter.add(bond_instance)
        self.system.part[0].bonds = ((bond_instance, 1),)

        # n+1 steps from min_dist to cut, then we remove the cut, because that
        # may break the bond due to rounding errors
        for dist in np.linspace(min_dist, cutoff, self.steps + 1)[:-1]:
            self.system.part[1].pos = self.system.part[
                0].pos + self.axis * dist
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["total"]
            E_ref = energy_func(dist)

            # Calculate forces
            f0_sim = np.copy(self.system.part[0].f)
            f1_sim = np.copy(self.system.part[1].f)
            f1_ref = self.axis * force_func(dist)

            # Check that energies match, ...
            self.assertAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force ...
            np.testing.assert_allclose(f0_sim, -f1_sim, 1E-12)
            # and has correct value.
            np.testing.assert_almost_equal(f1_sim, f1_ref)

            # Pressure
            # Isotropic pressure = 1/3 Trace Stress tensor
            # = 1/(3V) sum_i f_i r_i
            # where F is the force between the particles and r their distance
            p_expected = 1. / 3. * \
                np.dot(f1_sim, self.axis * dist) / self.system.volume()
            p_sim = self.system.analysis.pressure()["total"]
            self.assertAlmostEqual(p_sim, p_expected, delta=1E-12)

            # Pressure tensor
            # P_ij = 1/V F_i r_j
            p_tensor_expected = np.outer(
                f1_sim, self.axis * dist) / self.system.volume()
            p_tensor_sim = self.system.analysis.stress_tensor()["total"]
            np.testing.assert_allclose(
                p_tensor_sim, p_tensor_expected, atol=1E-12)
        if test_breakage:
            self.system.part[1].pos = self.system.part[0].pos \
                + self.axis * cutoff * (1.01)
            with self.assertRaisesRegex(Exception, "Encountered errors during integrate"):
                self.system.integrator.run(recalc_forces=True, steps=0)


if __name__ == '__main__':
    ut.main()
