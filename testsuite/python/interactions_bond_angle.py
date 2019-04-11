#
# Copyright (C) 2013-2018 The ESPResSo project
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
import espressomd
import numpy as np
import unittest as ut


class InteractionsAngleBondTest(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    box_l = 10.

    start_pos = np.random.rand(3) * box_l
    axis = np.random.rand(3)
    axis /= np.linalg.norm(axis)
    rel_pos = np.cross(np.random.rand(3), axis)
    rel_pos /= np.linalg.norm(rel_pos)
    N = 111  # number of angles to sample

    def setUp(self):
        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos + self.rel_pos, type=0)
        self.system.part.add(id=2, pos=self.start_pos + self.rel_pos, type=0)

        # Add a pair bond to make sure that doesn't cause trouble
        harmonic_bond = espressomd.interactions.HarmonicBond(k=0, r_0=0)
        self.system.bonded_inter.add(harmonic_bond)
        self.harmonic_bond = harmonic_bond
        self.system.part[1].add_bond((harmonic_bond, 0))

    def tearDown(self):
        self.system.part.clear()

    def rotate_vector(self, v, k, phi):
        """Rotates vector v around unit vector k by angle phi
        using Rodrigues' rotation formula.

        """
        vrot = v * np.cos(phi) + np.cross(k, v) * \
            np.sin(phi) + k * np.dot(k, v) * (1 - np.cos(phi))
        return vrot

    # Analytical expressions
    # Forces aren't divided by -sin(phi) to simplify the force magnitude check,
    # which only depends on the displacement vectors and potential gradient,
    # see "Bond angle potentials" in Doxygen
    def angle_harmonic_potential(self, phi, bend=1.0, phi0=np.pi):
        return 0.5 * bend * np.power(phi - phi0, 2)

    def angle_harmonic_force(self, phi, bend=1.0, phi0=np.pi):
        return bend * (phi - phi0)

    def angle_cosine_potential(self, phi, bend=1.0, phi0=np.pi):
        return bend * (1 - np.cos(phi - phi0))

    def angle_cosine_force(self, phi, bend=1.0, phi0=np.pi):
        return bend * np.sin(phi - phi0)

    def angle_cos_squared_potential(self, phi, bend=1.0, phi0=np.pi):
        return 0.5 * bend * (np.cos(phi) - np.cos(phi0))**2

    def angle_cos_squared_force(self, phi, bend=1.0, phi0=np.pi):
        return -bend * (np.cos(phi) - np.cos(phi0)) * np.sin(phi)

    def run_test(self, bond_instance, force_func, energy_func, phi0):
        self.system.bonded_inter.add(bond_instance)
        p0 = self.system.part[0]
        p1 = self.system.part[1]
        p2 = self.system.part[2]
        # Add bond angle
        p0.add_bond((bond_instance, 1, 2))
        # Add an extra (strength 0) pair bond, which should change nothing
        p0.add_bond((self.harmonic_bond, 1))

        d_phi = np.pi / self.N
        for i in range(1, self.N):  # avoid corner cases at phi = 0 and phi = pi
            p2.pos = self.start_pos + \
                self.rotate_vector(self.rel_pos, self.axis, i * d_phi)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = energy_func(i * d_phi)
            # Check that energies match
            np.testing.assert_almost_equal(E_sim, E_ref, decimal=4)

            f_ref = force_func(i * d_phi)
            phi_diff = i * d_phi - phi0
            for p, sign in [[p1, -1], [p2, +1]]:
                # Check that force is perpendicular
                self.assertAlmostEqual(
                    np.dot(p.f, self.system.distance_vec(p0, p)), 0,
                    delta=1E-12, msg="The force is not perpendicular")
                # Check that force has correct magnitude
                self.assertAlmostEqual(np.linalg.norm(p.f), np.abs(
                    f_ref) / self.system.distance(p0, p), delta=1E-11)
                # Check that force decreases the quantity abs(phi - phi0)
                if np.abs(phi_diff) > 1E-12:  # sign undefined for phi = phi0
                    force_axis = np.cross(self.system.distance_vec(p, p0), p.f)
                    self.assertEqual(
                        np.sign(phi_diff),
                        np.sign(sign * np.dot(force_axis, self.axis)),
                        msg="The force moves particles in the wrong direction")

            # Total force =0?
            np.testing.assert_allclose(
                np.sum(np.copy(self.system.part[0:3].f), 0), [0, 0, 0], atol=1E-12)

            # No pressure (isotropic compression preserves angles)
            self.assertAlmostEqual(
                self.system.analysis.pressure()["bonded"], 0, delta=1E-12)
            # Stress tensor trace=0 (isotropic compression preserves angles, 
            # consistency with pressure)
            self.assertAlmostEqual(
                np.trace(self.system.analysis.stress_tensor()["bonded"]), 0, delta=1E-12)

            # Pressure_ij =sum_p 1/V *F_i r_j
            # with r position of particle p and F force on particle p.
            # Then using F_p1=-F_p2 - F_p3 (no net force)
            # and r_p2 =r_p1 +r_{p1,p2} and r_p3 =r_p1 +r_{p1,p3}
            # P_ij =1/V (F_p2 r_{p1,p2} + F_p3 r_{p1,p3})
            p_tensor_expected = \
                np.outer(p1.f, self.system.distance_vec(p0, p1)) \
                + np.outer(p2.f, self.system.distance_vec(p0, p2))
            p_tensor_expected /= self.system.volume()
            np.testing.assert_allclose(
                self.system.analysis.stress_tensor()["bonded"],
                p_tensor_expected, atol=1E-12)

        # Remove bonds
        p0.delete_bond((bond_instance, 1, 2))
        p0.delete_bond((self.harmonic_bond, 1))
        p0.delete_bond((self.harmonic_bond, 2))

    def test_angle_harmonic(self):
        ah_bend = 1.
        ah_phi0 = 0.4327 * np.pi

        angle_harmonic = espressomd.interactions.AngleHarmonic(
            bend=ah_bend, phi0=ah_phi0)
        self.run_test(angle_harmonic,
                      lambda phi: self.angle_harmonic_force(
                      phi=phi, bend=ah_bend, phi0=ah_phi0),
                      lambda phi: self.angle_harmonic_potential(
                      phi=phi, bend=ah_bend, phi0=ah_phi0),
                      ah_phi0)

    def test_angle_cosine(self):
        ac_bend = 1
        ac_phi0 = 1

        angle_cosine = espressomd.interactions.AngleCosine(
            bend=ac_bend, phi0=ac_phi0)
        self.run_test(angle_cosine,
                      lambda phi: self.angle_cosine_force(
                      phi=phi, bend=ac_bend, phi0=ac_phi0),
                      lambda phi: self.angle_cosine_potential(
                      phi=phi, bend=ac_bend, phi0=ac_phi0),
                      ac_phi0)

    def test_angle_cos_squared(self):
        acs_bend = 1
        acs_phi0 = 1

        angle_cos_squared = espressomd.interactions.AngleCossquare(
            bend=acs_bend, phi0=acs_phi0)
        self.run_test(angle_cos_squared,
                      lambda phi: self.angle_cos_squared_force(
                      phi=phi, bend=acs_bend, phi0=acs_phi0),
                      lambda phi: self.angle_cos_squared_potential(
                      phi=phi, bend=acs_bend, phi0=acs_phi0),
                      acs_phi0)

    @ut.skipIf(not espressomd.has_features("TABULATED"),
               "Skipped because feature is disabled")
    def test_angle_tabulated(self):
        """Check that we can reproduce the three other potentials."""
        at_bend = 1
        at_phi0 = 1
        phi = np.linspace(0, np.pi, num=self.N + 1, endpoint=True)
        for fun_pot, fun_force in [
                (self.angle_cosine_potential, self.angle_cosine_force),
                (self.angle_harmonic_potential, self.angle_harmonic_force),
                (self.angle_cos_squared_potential, self.angle_cos_squared_force)]:
            angle_tabulated = espressomd.interactions.Tabulated(
                type='angle', min=0, max=np.pi,
                energy=fun_pot(phi=phi, bend=at_bend, phi0=at_phi0),
                force=fun_force(phi=phi, bend=at_bend, phi0=at_phi0))
            self.run_test(angle_tabulated,
                          lambda phi: fun_force(
                          phi=phi, bend=at_bend, phi0=at_phi0),
                          lambda phi: fun_pot(
                          phi=phi, bend=at_bend, phi0=at_phi0),
                          at_phi0)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
