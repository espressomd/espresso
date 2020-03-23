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
import numpy as np

import espressomd

# Dihedral interaction needs more rigorous tests.
# The geometry checked here is rather simple and special.
# I also found that as the dihedral angle approaches to 0, the simulation
# values deviate from the analytic values by roughly 10%.


def rotate_vector(v, k, phi):
    """Rotates vector v around unit vector k by angle phi.
    Uses Rodrigues' rotation formula."""
    vrot = v * np.cos(phi) + np.cross(k, v) * \
        np.sin(phi) + k * np.dot(k, v) * (1.0 - np.cos(phi))
    return vrot


def dihedral_potential(k, phi, n, phase):
    if phi == -1:
        return 0
    else:
        return k * (1 - np.cos(n * phi - phase))


def dihedral_force(k, n, phase, p1, p2, p3, p4):
    v12 = p2 - p1
    v23 = p3 - p2
    v34 = p4 - p3

    v12Xv23 = np.cross(v12, v23)
    l_v12Xv23 = np.linalg.norm(v12Xv23)
    v23Xv34 = np.cross(v23, v34)
    l_v23Xv34 = np.linalg.norm(v23Xv34)
    # if dihedral angle is not defined, no forces
    if l_v12Xv23 <= 1e-8 or l_v23Xv34 <= 1e-8:
        return 0, 0, 0
    else:
        cosphi = np.abs(np.dot(v12Xv23, v23Xv34)) / (l_v12Xv23 * l_v23Xv34)
        phi = np.arccos(cosphi)
        f1 = (v23Xv34 - cosphi * v12Xv23) / l_v12Xv23
        f4 = (v12Xv23 - cosphi * v23Xv34) / l_v23Xv34

        v23Xf1 = np.cross(v23, f1)
        v23Xf4 = np.cross(v23, f4)
        v34Xf4 = np.cross(v34, f4)
        v12Xf1 = np.cross(v12, f1)

        coeff = -k * n * np.sin(n * phi - phase) / np.sin(phi)

        force1 = coeff * v23Xf1
        force2 = coeff * (v34Xf4 - v12Xf1 - v23Xf1)
        force3 = coeff * (v12Xf1 - v23Xf4 - v34Xf4)
        return force1, force2, force3


class InteractionsBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(seed=42)

    box_l = 10.

    start_pos = [5., 5., 5.]
    axis = np.array([1., 0., 0.])
    axis /= np.linalg.norm(axis)
    rel_pos_1 = np.array([0., 1., 0.])
    rel_pos_2 = np.array([0., 0., 1.])

    def setUp(self):

        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos, type=0)
        self.system.part.add(id=2, pos=self.start_pos, type=0)
        self.system.part.add(id=3, pos=self.start_pos, type=0)

    def tearDown(self):
        self.system.part.clear()

    # Analytical Expression
    def dihedral_angle(self, p1, p2, p3, p4):
        """
        Calculate the dihedral angle phi based on particles' position p1, p2, p3, p4.
        """
        v12 = p2 - p1
        v23 = p3 - p2
        v34 = p4 - p3

        v12Xv23 = np.cross(v12, v23)
        l_v12Xv23 = np.linalg.norm(v12Xv23)
        v23Xv34 = np.cross(v23, v34)
        l_v23Xv34 = np.linalg.norm(v23Xv34)

        # if dihedral angle is not defined, phi := -1.
        if l_v12Xv23 <= 1e-8 or l_v23Xv34 <= 1e-8:
            return -1
        else:
            cosphi = np.abs(np.dot(v12Xv23, v23Xv34)) / (
                l_v12Xv23 * l_v23Xv34)
            return np.arccos(cosphi)

    # Test Dihedral Angle
    def test_dihedral(self):
        dh_k = 1
        dh_phase = np.pi / 6
        dh_n = 1

        dh = espressomd.interactions.Dihedral(
            bend=dh_k, mult=dh_n, phase=dh_phase)
        self.system.bonded_inter.add(dh)
        self.system.part[1].add_bond((dh, 0, 2, 3))
        self.system.part[2].pos = self.system.part[1].pos + [1, 0, 0]

        N = 111
        d_phi = np.pi / (N * 4)
        for i in range(N):
            self.system.part[0].pos = self.system.part[1].pos + \
                rotate_vector(self.rel_pos_1, self.axis, i * d_phi)
            self.system.part[3].pos = self.system.part[2].pos + \
                rotate_vector(self.rel_pos_2, self.axis, -i * d_phi)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            phi = self.dihedral_angle(self.system.part[0].pos,
                                      self.system.part[1].pos,
                                      self.system.part[2].pos,
                                      self.system.part[3].pos)
            E_ref = dihedral_potential(dh_k, phi, dh_n, dh_phase)

            # Calculate forces
            f2_sim = self.system.part[1].f
            _, f2_ref, _ = dihedral_force(dh_k, dh_n, dh_phase,
                                          self.system.part[0].pos,
                                          self.system.part[1].pos,
                                          self.system.part[2].pos,
                                          self.system.part[3].pos)

            # Check that energies match, ...
            np.testing.assert_almost_equal(E_sim, E_ref)
            # and has correct value.
            f2_sim_copy = np.copy(f2_sim)
            np.testing.assert_almost_equal(f2_sim_copy, f2_ref)

    # Test Tabulated Dihedral Angle
    def test_tabulated_dihedral(self):
        N = 111
        d_phi = 2 * np.pi / N
        # tabulated values for the range [0, 2*pi]
        tab_energy = [np.cos(i * d_phi) for i in range(N + 1)]
        tab_force = [np.cos(i * d_phi) for i in range(N + 1)]

        dihedral_tabulated = espressomd.interactions.TabulatedDihedral(
            energy=tab_energy, force=tab_force)
        self.system.bonded_inter.add(dihedral_tabulated)
        self.system.part[1].add_bond((dihedral_tabulated, 0, 2, 3))
        self.system.part[2].pos = self.system.part[1].pos + [1, 0, 0]

        # check stored parameters
        interaction_id = len(self.system.bonded_inter) - 1
        tabulated = self.system.bonded_inter[interaction_id]
        np.testing.assert_allclose(tabulated.params['force'], tab_force)
        np.testing.assert_allclose(tabulated.params['energy'], tab_energy)
        np.testing.assert_almost_equal(tabulated.params['min'], 0.)
        np.testing.assert_almost_equal(tabulated.params['max'], 2 * np.pi)

        # measure at half the angular resolution to observe interpolation
        for i in range(2 * N - 1):
            # increase dihedral angle by d_phi (phi ~ 0 at i = 0)
            self.system.part[0].pos = self.system.part[1].pos + \
                rotate_vector(self.rel_pos_1, self.axis, -i * d_phi / 4)
            self.system.part[3].pos = self.system.part[2].pos + \
                rotate_vector(self.rel_pos_1, self.axis, i * d_phi / 4)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]

            # Get tabulated values
            j = i // 2
            if i % 2 == 0:
                E_ref = tab_energy[j]
            else:
                E_ref = (tab_energy[j] + tab_energy[j + 1]) / 2.0

            # Check that energies match, ...
            np.testing.assert_almost_equal(E_sim, E_ref)


if __name__ == '__main__':
    ut.main()
