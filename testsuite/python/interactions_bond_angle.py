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


@ut.skipIf(not espressomd.has_features(["BOND_ANGLE"]),
           "Features not available, skipping test!")
class InteractionsNonBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    box_l = 10.

    start_pos = np.random.rand(3) * box_l
    axis = np.random.rand(3)
    axis /= np.linalg.norm(axis)
    rel_pos = np.cross(np.random.rand(3), axis)
    rel_pos /= np.linalg.norm(rel_pos)

    def setUp(self):
        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos + self.rel_pos, type=0)
        self.system.part.add(id=2, pos=self.start_pos + self.rel_pos, type=0)

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
    def angle_harmonic_potential(self, phi, bend=1.0, phi0=np.pi):
        return 0.5 * bend * np.power(phi - phi0, 2)

    def angle_cosine_potential(self, phi, bend=1.0, phi0=np.pi):
        return bend * (1 - np.cos(phi - phi0))

    def angle_cos_squared_potential(self, phi, bend=1.0, phi0=np.pi):
        return 0.5 * bend * (np.cos(phi) - np.cos(phi0))**2

    def test_angle_harmonic(self):
        ah_bend = 1.
        ah_phi0 = 0.4327 * np.pi

        angle_harmonic = espressomd.interactions.AngleHarmonic(
            bend=ah_bend, phi0=ah_phi0)
        self.system.bonded_inter.add(angle_harmonic)
        self.system.part[0].add_bond((angle_harmonic, 1, 2))

        N = 111
        d_phi = np.pi / N
        for i in range(N):
            self.system.part[2].pos = self.start_pos + \
                self.rotate_vector(self.rel_pos, self.axis, i * d_phi)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = self.angle_harmonic_potential(
                phi=i * d_phi, bend=ah_bend, phi0=ah_phi0)

            # Calculate forces
            f0_sim = np.copy(self.system.part[0].f)
            f1_sim = np.copy(self.system.part[1].f)
            f2_sim = np.copy(self.system.part[2].f)

            # Check that energies match
            np.testing.assert_almost_equal(E_sim, E_ref, decimal=4)
            # Force equals minus the counter-force
            np.testing.assert_allclose(f0_sim + f2_sim, -f1_sim)

    # Test Angle Cosine Potential
    def test_angle_cosine(self):
        ac_bend = 1
        ac_phi0 = 1

        angle_cosine = espressomd.interactions.AngleCosine(
            bend=ac_bend, phi0=ac_phi0)
        self.system.bonded_inter.add(angle_cosine)
        self.system.part[0].add_bond((angle_cosine, 1, 2))

        N = 111
        d_phi = np.pi / N
        for i in range(N):
            self.system.part[2].pos = self.start_pos + \
                self.rotate_vector(self.rel_pos, self.axis, i * d_phi)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()['bonded']
            E_ref = self.angle_cosine_potential(
                phi=i * d_phi, bend=ac_bend, phi0=ac_phi0)

            # Calculate forces
            f0_sim = np.copy(self.system.part[0].f)
            f1_sim = np.copy(self.system.part[1].f)
            f2_sim = np.copy(self.system.part[2].f)

            # Check that energies match
            np.testing.assert_almost_equal(E_sim, E_ref, decimal=4)
            # Force equals minus the counter-force
            np.testing.assert_allclose(f0_sim + f2_sim, -f1_sim)

    def test_angle_cos_squared(self):
        acs_bend = 1
        acs_phi0 = 1

        angle_cos_squared = espressomd.interactions.AngleCossquare(
            bend=acs_bend, phi0=acs_phi0)
        self.system.bonded_inter.add(angle_cos_squared)
        self.system.part[0].add_bond((angle_cos_squared, 1, 2))

        N = 111
        d_phi = np.pi / N
        for i in range(N):
            self.system.part[2].pos = self.start_pos + \
                self.rotate_vector(self.rel_pos, self.axis, i * d_phi)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()['bonded']
            E_ref = self.angle_cos_squared_potential(
                phi=i * d_phi, bend=acs_bend, phi0=acs_phi0)

            # Calculate forces
            f0_sim = np.copy(self.system.part[0].f)
            f1_sim = np.copy(self.system.part[1].f)
            f2_sim = np.copy(self.system.part[2].f)

            # Check that energies match
            np.testing.assert_almost_equal(E_sim, E_ref)
            # Force equals minus the counter-force
            np.testing.assert_allclose(f0_sim + f2_sim, -f1_sim)

    # Test Tabulated Harmonic Angle
    @ut.skipIf(not espressomd.has_features(["TABULATED"]),
               "TABULATED feature is not available, skipping tabulated test.")
    def test_angle_tabulated(self):
        ah_bend = 1.
        ah_phi0 = 0.4327 * np.pi
        N = 111
        d_phi = np.pi / N

        tab_energy = [self.angle_harmonic_potential(
            phi=i * d_phi, bend=ah_bend, phi0=ah_phi0) for i in range(N + 1)]
        tab_force = [ah_bend * abs(i * d_phi - ah_phi0) for i in range(N + 1)]
        angle_tabulated = espressomd.interactions.Tabulated(
            type='angle', energy=tab_energy, force=tab_force, min=0., max=np.pi)
        self.system.bonded_inter.add(angle_tabulated)
        self.system.part[0].add_bond((angle_tabulated, 1, 2))

        # check stored parameters
        interaction_id = len(self.system.bonded_inter) - 1
        tabulated = self.system.bonded_inter[interaction_id]
        np.testing.assert_allclose(tabulated.params['force'], tab_force)
        np.testing.assert_allclose(tabulated.params['energy'], tab_energy)
        np.testing.assert_almost_equal(tabulated.params['min'], 0.)
        np.testing.assert_almost_equal(tabulated.params['max'], np.pi)

        # measure at half the angular resolution to observe interpolation
        for i in range(2 * N - 1):
            self.system.part[2].pos = self.start_pos + \
                self.rotate_vector(self.rel_pos, self.axis, i * d_phi / 2.0)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]

            # Calculate forces
            f0_sim = np.copy(self.system.part[0].f)
            f1_sim = np.copy(self.system.part[1].f)
            f2_sim = np.copy(self.system.part[2].f)

            # Reference value: traverse the tabulated array in reverse using a
            # negative index (ESPResSo uses the external angle convention)
            j = -(i // 2) - 1
            if i % 2 == 0:
                E_ref = tab_energy[j]
                f0_ref = tab_force[j]
            else:
                E_ref = (tab_energy[j] + tab_energy[j - 1]) / 2.0
                f0_ref = (tab_force[j] + tab_force[j - 1]) / 2.0

            # Check that energies match, ...
            np.testing.assert_almost_equal(E_sim, E_ref, decimal=4)
            # Force equals minus the counter-force
            np.testing.assert_allclose(f0_sim + f2_sim, -f1_sim)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
