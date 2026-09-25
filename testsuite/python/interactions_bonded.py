#
# Copyright (C) 2013-2026 The ESPResSo project
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
import espressomd.interactions

#
# Analytical expressions for interactions
#


def harmonic_potential(scalar_r, k, r_0):
    return 0.5 * k * (scalar_r - r_0)**2


def harmonic_force(scalar_r, k, r_0):
    return -k * (scalar_r - r_0)


def fene_potential(scalar_r, k, d_r_max, r_0):
    return -0.5 * k * d_r_max**2 * np.log(1 - ((scalar_r - r_0) / d_r_max)**2)


def fene_force(scalar_r, k, d_r_max, r_0):
    return k * (scalar_r - r_0) * d_r_max**2 / \
        ((scalar_r - r_0)**2 - d_r_max**2)


def coulomb_potential(scalar_r, k, q1, q2):
    return k * q1 * q2 / scalar_r


def coulomb_force(scalar_r, k, q1, q2):
    return k * q1 * q2 / scalar_r**2


def quartic_force(k0, k1, r, r_cut, scalar_r):
    force = 0.
    if scalar_r <= r_cut:
        force = - k0 * (scalar_r - r) - k1 * (scalar_r - r)**3
    return force


def quartic_potential(k0, k1, r, r_cut, scalar_r):
    energy = 0.
    if scalar_r <= r_cut:
        energy = 0.5 * k0 * (scalar_r - r)**2 + 0.25 * k1 * (scalar_r - r)**4
    return energy


#
# Test case
#

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

        self.system.part.add(pos=self.start_pos, type=0)
        self.system.part.add(pos=self.start_pos, type=0)

    def tearDown(self):
        if espressomd.has_features(["ELECTROSTATICS"]):
            self.system.electrostatics.clear()
        if espressomd.has_features(['LENNARD_JONES']):
            self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()
        self.system.part.clear()
        self.system.bonded_inter.clear()

    # Test Harmonic Bond
    def test_harmonic(self):
        k = 5.
        r_0 = 1.5
        r_cut = 3.355

        hb = espressomd.interactions.HarmonicBond(k=k, r_0=r_0, r_cut=r_cut)
        self.run_test(hb,
                      lambda r: harmonic_force(r, k, r_0),
                      lambda r: harmonic_potential(r, k, r_0),
                      0.01, r_cut, True, test_same_pos_exception=True)

    # Test Fene Bond
    def test_fene(self):
        k = 23.15
        d_r_max = 3.355
        r_0 = 1.1

        fene = espressomd.interactions.FeneBond(k=k, d_r_max=d_r_max, r_0=r_0)
        self.run_test(fene,
                      lambda r: fene_force(r, k, d_r_max, r_0),
                      lambda r: fene_potential(r, k, d_r_max, r_0=r_0),
                      0.01, r_0 + d_r_max, True, test_same_pos_exception=True)

    def test_virtual_bond(self):
        # add sentinel harmonic bond, otherwise short-range loop is skipped
        hb = espressomd.interactions.HarmonicBond(k=1., r_0=0.1, r_cut=0.5)
        vb = espressomd.interactions.Virtual()
        self.system.bonded_inter.add(hb)
        self.system.bonded_inter.add(vb)
        p1, p2 = self.system.part.all()
        p1.add_bond((vb, p2))

        self.system.integrator.run(steps=0, recalc_forces=True)
        self.assertEqual(self.system.analysis.energy()["total"], 0.)
        np.testing.assert_allclose(np.copy(p1.f), 0., atol=1e-12, rtol=0)
        np.testing.assert_allclose(np.copy(p2.f), 0., atol=1e-12, rtol=0)

    @utx.skipIfMissingFeatures(["BOND_CONSTRAINT"])
    def test_rigid_bond(self):
        rb = espressomd.interactions.RigidBond(r=1.0, ptol=0.1, vtol=0.1)
        self.system.bonded_inter.add(rb)
        p1, p2 = self.system.part.all()
        p2.pos = p1.pos + np.array([1.0, 0., 0.])
        p1.add_bond((rb, p2))

        self.system.integrator.run(steps=0, recalc_forces=True)
        self.assertEqual(self.system.analysis.energy()["total"], 0.)
        np.testing.assert_allclose(np.copy(p1.f), 0., atol=1e-12, rtol=0)
        np.testing.assert_allclose(np.copy(p2.f), 0., atol=1e-12, rtol=0)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_coulomb(self):
        k = 1.
        q1 = 1.
        q2 = -1.
        p1, p2 = self.system.part.all()
        p1.q = q1
        p2.q = q2
        self.run_test(
            espressomd.interactions.BondedCoulomb(prefactor=k),
            lambda r: coulomb_force(r, k, q1, q2),
            lambda r: coulomb_potential(r, k, q1, q2),
            0.01, self.system.box_l[0] / 3)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_coulomb_ghost_partner(self):
        """Regression test: a bonded-Coulomb pair whose partner is reached
        through the deduped ghost tail. With a 2-rank ``[2, 1, 1]`` node grid
        the box is split at ``x = box_l / 2``; placing the two partners on
        opposite sides of that plane forces the bond partner to be resolved as
        a ghost, i.e. at a pack index in the tail where ``row(idx) != idx``.
        The bonded-Coulomb energy/pressure/force kernels must translate the
        pack index to the store row before reading the charge; reading the
        charge by raw pack index yields the wrong store row (and hence the
        wrong ``q1 * q2`` product) for the ghost partner.
        """
        system = self.system
        system.part.clear()
        k = 1.
        q1 = 1.
        q2 = -1.
        # Straddle the x = box_l / 2 rank boundary of the [2, 1, 1] node grid.
        half = system.box_l[0] / 2.
        dist = 1.0
        p1 = system.part.add(pos=[half - dist / 2., 5., 5.], q=q1, type=0)
        p2 = system.part.add(pos=[half + dist / 2., 5., 5.], q=q2, type=0)
        axis = np.array([1., 0., 0.])

        bond = espressomd.interactions.BondedCoulomb(prefactor=k)
        system.bonded_inter.add(bond)
        p1.bonds = ((bond, p2),)

        system.integrator.run(recalc_forces=True, steps=0)

        # Energy: q1 * q2 = -1 must survive the store-row translation.
        E_sim = system.analysis.energy()["total"]
        E_ref = coulomb_potential(dist, k, q1, q2)
        self.assertAlmostEqual(E_sim, E_ref, places=7)

        # Forces: equal and opposite, attractive along the separation axis.
        f0_sim = np.copy(p1.f)
        f1_sim = np.copy(p2.f)
        f1_ref = axis * coulomb_force(dist, k, q1, q2)
        np.testing.assert_allclose(f0_sim, -f1_sim, atol=1e-12)
        np.testing.assert_allclose(f1_sim, f1_ref, atol=1e-12)

        # Isotropic pressure = 1/(3V) sum_i f_i . r_i.
        p_expected = np.dot(f1_sim, axis * dist) / (3. * system.volume())
        p_sim = system.analysis.pressure()["total"]
        self.assertAlmostEqual(p_sim, p_expected, delta=1e-12)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_coulomb_sr(self):
        # with negated actual charges and only short range int: cancels out all
        # interactions
        q1 = 1.2
        q2 = -q1
        p1, p2 = self.system.part.all()
        p1.q = q1
        p2.q = q2
        r_cut = 2

        sr_solver = espressomd.electrostatics.DH(
            prefactor=2, kappa=0.8, r_cut=r_cut)
        self.system.electrostatics.solver = sr_solver
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
        k0 = 2.
        k1 = 5.
        q_r = 0.5
        r_cut = self.system.box_l[0] / 3.

        quartic = espressomd.interactions.QuarticBond(k0=k0,
                                                      k1=k1,
                                                      r=q_r,
                                                      r_cut=r_cut)

        self.run_test(quartic,
                      lambda r: quartic_force(k0, k1, q_r, r_cut, r),
                      lambda r: quartic_potential(k0, k1, q_r, r_cut, r),
                      0.01, r_cut, True, test_same_pos_exception=True)

    def run_test(self, bond_instance, force_func, energy_func, min_dist,
                 cutoff, test_breakage=False, test_same_pos_exception=False):
        self.system.bonded_inter.add(bond_instance)
        p1, p2 = self.system.part.all()
        p1.bonds = ((bond_instance, p2),)

        # n+1 steps from min_dist to cut, then we remove the cut, because that
        # may break the bond due to rounding errors
        for dist in np.linspace(min_dist, cutoff, self.steps + 1)[:-1]:
            p2.pos = p1.pos + self.axis * dist
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["total"]
            E_ref = energy_func(dist)

            # Calculate forces
            f0_sim = np.copy(p1.f)
            f1_sim = np.copy(p2.f)
            f1_ref = self.axis * force_func(dist)

            # Check that energies match, ...
            self.assertAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force ...
            np.testing.assert_allclose(f0_sim, -f1_sim, 1E-12)
            # and has correct value.
            np.testing.assert_almost_equal(f1_sim, f1_ref)

            # Pressure
            # Isotropic pressure = 1/3 trace(pressure tensor)
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
            p_tensor_sim = self.system.analysis.pressure_tensor()["total"]
            np.testing.assert_allclose(
                p_tensor_sim, p_tensor_expected, atol=1E-12)
        if test_breakage:
            p2.pos = p1.pos + self.axis * cutoff * 1.01
            with self.assertRaisesRegex(Exception, r"while calling method integrate\(\)"):
                self.system.integrator.run(recalc_forces=True, steps=0)
            with self.assertRaisesRegex(Exception, r"while calling method calculate_energy\(\)"):
                self.system.analysis.energy()["total"]
        if test_same_pos_exception:
            p2.pos = p1.pos
            with self.assertRaises(Exception):
                self.system.integrator.run(recalc_forces=True, steps=0)

    def test_fene_bond_breakage(self):
        system = self.system
        system.part.clear()
        p1 = system.part.add(pos=system.box_l / 2. + [0., 0., -2.], id=1)
        p2 = system.part.add(pos=system.box_l / 2. + [0., 0., -1.], id=2)
        p3 = system.part.add(pos=system.box_l / 2. + [0., 0., +1.], id=3)
        p4 = system.part.add(pos=system.box_l / 2. + [0., 0., +2.], id=4)

        # Allow bond length in the interval (1.5, 2.5)
        epsilon = 1e-6
        fene = espressomd.interactions.FeneBond(
            k=1., d_r_max=0.5 - epsilon, r_0=2.)
        system.bonded_inter.add(fene)

        # Expect no bond breakage
        p2.bonds = (fene, p3)  # dist 2.0
        system.integrator.run(0, recalc_forces=True)
        p2.delete_all_bonds()

        # Expect bond breakage
        error_msg = "ERROR: bond broken between particles {}, {}"

        with self.subTest(msg="FENE bond breaks on elongation"):
            p1.bonds = (fene, p4)  # dist 4.0
            with self.assertRaisesRegex(Exception, error_msg.format(1, 4)):
                system.integrator.run(0, recalc_forces=True)
            with self.assertRaisesRegex(Exception, error_msg.format(1, 4)):
                self.assertEqual(self.system.analysis.energy()["total"], 0.)
            p1.delete_all_bonds()

            p2.bonds = (fene, p4)  # dist 3.0
            with self.assertRaisesRegex(Exception, error_msg.format(2, 4)):
                system.integrator.run(0, recalc_forces=True)
            with self.assertRaisesRegex(Exception, error_msg.format(2, 4)):
                self.assertEqual(self.system.analysis.energy()["total"], 0.)
            p2.delete_all_bonds()

        with self.subTest(msg="FENE bond breaks on compression"):
            p1.bonds = (fene, p2)  # dist 1.0
            # Expect bond breakage
            with self.assertRaisesRegex(Exception, error_msg.format(1, 2)):
                system.integrator.run(0, recalc_forces=True)
            with self.assertRaisesRegex(Exception, error_msg.format(1, 2)):
                self.assertEqual(self.system.analysis.energy()["total"], 0.)
            p1.delete_all_bonds()

            p3.bonds = (fene, p4)  # dist 1.0
            # Expect bond breakage
            with self.assertRaisesRegex(Exception, error_msg.format(3, 4)):
                system.integrator.run(0, recalc_forces=True)
            with self.assertRaisesRegex(Exception, error_msg.format(3, 4)):
                self.assertEqual(self.system.analysis.energy()["total"], 0.)
            p3.delete_all_bonds()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] < 2,
               "only runs for 2 or more MPI ranks")
    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_bond_resolution_error(self):
        system = self.system
        system.part.clear()

        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=1.0, shift="auto")

        hb = espressomd.interactions.HarmonicBond(k=1.0, r_0=1.0, r_cut=0.0)
        system.bonded_inter.add(hb)

        system.part.add(id=0, pos=[2.5, 5.0, 5.0])
        system.part.add(id=1, pos=[7.5, 5.0, 5.0])
        system.part.by_id(0).add_bond((hb, 1))

        # Expect bond resolution error
        error_msg = "ERROR: bond partner not found on local node, could only find: {}"
        with self.assertRaisesRegex(Exception, error_msg.format(1)):
            system.integrator.run(0)


if __name__ == '__main__':
    ut.main()
