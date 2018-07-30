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
import espressomd
import numpy as np
import unittest as ut


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
        """Rotates vector v around unit vector k by angle phi.
        Uses Rodrigues' rotation formula."""
        vrot = v * np.cos(phi) + np.cross(k, v) * \
            np.sin(phi) + k * np.dot(k, v) * (1 - np.cos(phi))
        return vrot

    # Required, since assertAlmostEqual does NOT check significant places
    def assertFractionAlmostEqual(self, a, b, places=10):
        if abs(b) < 1E-8:
            self.assertAlmostEqual(a, b)
        else:
            self.assertAlmostEqual(a / b, 1., places=4)

    def assertItemsFractionAlmostEqual(self, a, b):
        for i, ai in enumerate(a):
            self.assertFractionAlmostEqual(ai, b[i])

    # Analytical Expression
    # for Angle Harmonic Potential
    def angle_harmonic_potential(self, phi, bend=1.0, phi0=np.pi):
        return 0.5 * bend * np.power(phi - phi0, 2)
    # for Angle Cosine Potential
    def angle_cosine_potential(self, phi, bend=1.0, phi0 = np.pi):
        return bend*(1-np.cos(phi - phi0))
    # for Angle Cosine Squared Potential
    def angle_cos_squared_potential(self, phi, bend = 1.0, phi0 = np.pi):
        return 0.5*bend*(np.cos(phi) - np.cos(phi0))**2

    # Test Angle Harmonic Potential
    @ut.skipIf(not espressomd.has_features(["BOND_ANGLE"]),
               "Features not available, skipping test!")
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

            # Check that energies match
            self.assertFractionAlmostEqual(E_sim, E_ref)

    # Test Angle Cosine Potential
    def test_angle_cosine(self):

        ac_bend = 1
        ac_phi0 = 1

        angle_cosine = espressomd.interactions.AngleCosine(
            bend = ac_bend, phi0 = ac_phi0)
        self.system.bonded_inter.add(angle_cosine)
        self.system.part[0].add_bond((angle_cosine, 1, 2))

        N = 111
        d_phi = np.pi / N
        for i in range(N):
            self.system.part[2].pos = self.start_pos + \
                    self.rotate_vector(self.rel_pos, self.axis, i*d_phi)
            self.system.integrator.run(recalc_forces = True, steps = 0)

            # Calculate energies
            E_sim = self.system.analysis.energy()['bonded']
            E_ref = self.angle_cosine_potential(
                phi = i*d_phi, bend = ac_bend, phi0 = ac_phi0)

            # Check that energies match
            self.assertFractionAlmostEqual(E_sim, E_ref)

    # Test Angle Cosine Squared Potential
    def test_angle_cos_squared(self):

        acs_bend = 1
        acs_phi0 = 1

        angle_cos_squared = espressomd.interactions.AngleCossquare(
            bend = acs_bend, phi0 = acs_phi0)
        self.system.bonded_inter.add(angle_cos_squared)
        self.system.part[0].add_bond((angle_cos_squared, 1, 2))

        N = 111
        d_phi = np.pi / N
        for i in range(N):
            self.system.part[2].pos = self.start_pos + \
                    self.rotate_vector(self.rel_pos, self.axis, i*d_phi)
            self.system.integrator.run(recalc_forces = True, steps = 0)

            # Calculate energies
            E_sim = self.system.analysis.energy()['bonded']
            E_ref = self.angle_cos_squared_potential(
                phi = i*d_phi, bend = acs_bend, phi0 = acs_phi0)

            # Check that energies match
            self.assertFractionAlmostEqual(E_sim, E_ref)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
