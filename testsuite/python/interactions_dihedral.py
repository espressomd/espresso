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


def rotate_vector(v, k, phi):
    """Rotates vector v around unit vector k by angle phi.
    Uses Rodrigues' rotation formula."""
    vrot = np.array(v) * np.cos(phi) + np.cross(k, v) * \
        np.sin(phi) + np.array(k) * np.dot(k, v) * (1.0 - np.cos(phi))
    return vrot


def dihedral_potential_and_forces(k, n, phase, p1, p2, p3, p4):
    """
    Calculate the potential and forces for a dihedral angle.
    """
    v12 = p2 - p1
    v23 = p3 - p2
    v34 = p4 - p3

    v12Xv23 = np.cross(v12, v23)
    l_v12Xv23 = np.linalg.norm(v12Xv23)
    v23Xv34 = np.cross(v23, v34)
    l_v23Xv34 = np.linalg.norm(v23Xv34)

    phi = np.arctan2(np.dot(v23, np.cross(v12Xv23, v23Xv34)),
                     np.dot(v23, v23) * np.dot(v12Xv23, v23Xv34))

    f1 = (v23Xv34 - np.cos(phi) * v12Xv23) / l_v12Xv23
    f4 = (v12Xv23 - np.cos(phi) * v23Xv34) / l_v23Xv34

    v23Xf1 = np.cross(v23, f1)
    v23Xf4 = np.cross(v23, f4)
    v34Xf4 = np.cross(v34, f4)
    v12Xf1 = np.cross(v12, f1)

    # handle singularity near TINY_SIN_VALUE
    if np.abs(np.sin(phi)) < 1e-10:
        coeff = -k * n**2 * np.cos(n * phi - phase) / np.cos(phi)
    else:
        coeff = -k * n * np.sin(n * phi - phase) / np.sin(phi)

    force1 = coeff * v23Xf1
    force2 = coeff * (v34Xf4 - v12Xf1 - v23Xf1)
    force3 = coeff * (v12Xf1 - v23Xf4 - v34Xf4)
    force4 = coeff * v23Xf4
    potential = k * (1 - np.cos(n * phi - phase))
    return (potential, (force1, force2, force3, force4))


class InteractionsBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.4
    system.time_step = 0.1
    np.random.seed(seed=42)

    def tearDown(self):
        self.system.part.clear()

    def check_values(self, E_ref, forces_ref, tol=1e-12):
        E_sim = self.system.analysis.energy()["bonded"]
        np.testing.assert_allclose(E_sim, E_ref, atol=tol)
        if forces_ref:
            f0, f1, f2, f3 = self.system.part.all().f
            f0_ref, f1_ref, f2_ref, f3_ref = forces_ref
            np.testing.assert_allclose(np.copy(f0), f0_ref, atol=tol)
            np.testing.assert_allclose(np.copy(f1), f1_ref, atol=tol)
            np.testing.assert_allclose(np.copy(f2), f2_ref, atol=tol)
            np.testing.assert_allclose(np.copy(f3), f3_ref, atol=tol)

    # Test Dihedral Angle
    def test_dihedral(self):
        axis = np.array([1., 0., 0.])
        dh_k = 2.
        N = 100  # even number to get singularities at phi=0 and phi=pi
        d_phi = 2 * np.pi / N
        for dh_n, dh_phi0_div in [(2, 3), (3, 6)]:
            with self.subTest(multiplicity=dh_n, phi_0=f"pi / {dh_phi0_div}"):
                dh_phi0 = np.pi / dh_phi0_div
                dihedral = espressomd.interactions.Dihedral(
                    bend=dh_k, mult=dh_n, phase=dh_phi0)
                self.system.bonded_inter.add(dihedral)
                self.system.part.clear()
                p0, p1, p2, p3 = self.system.part.add(pos=4 * [(0., 0., 0.)])
                p1.add_bond((dihedral, p0, p2, p3))
                p1.pos = [5., 5., 5.]
                p2.pos = p1.pos + [1., 0., 0.]
                p0.pos = p1.pos + [0., 1., 0.]

                for i in range(N):
                    phi = i * d_phi
                    p3.pos = p2.pos + rotate_vector([0., 1., 0.], axis, phi)
                    self.system.integrator.run(recalc_forces=True, steps=0)

                    # Calculate expected forces and energies
                    E_ref, forces_ref = dihedral_potential_and_forces(
                        dh_k, dh_n, dh_phi0, p0.pos, p1.pos, p2.pos, p3.pos)

                    self.check_values(E_ref, forces_ref)

    # Test Tabulated Dihedral Angle
    @utx.skipIfMissingFeatures(["TABULATED"])
    def test_tabulated_dihedral(self):
        axis = np.array([1., 0., 0.])
        dh_k = 2.
        N = 100  # even number to get singularities at phi=0 and phi=pi
        d_phi = 2 * np.pi / N
        for dh_n, dh_phi0_div in [(2, 3), (3, 6)]:
            with self.subTest(multiplicity=dh_n, phi_0=f"pi / {dh_phi0_div}"):
                dh_phi0 = np.pi / dh_phi0_div
                # tabulated values for the range [0, 2*pi]
                phi = d_phi * np.arange(N + 1)
                tab_energy = dh_k * (1. - np.cos(dh_n * phi - dh_phi0))
                div = np.sin(phi)
                div[0] = div[N // 2] = div[N] = 1.
                tab_force = -dh_k * dh_n * np.sin(dh_n * phi - dh_phi0) / div
                tab_force[0] = tab_force[N // 2] = tab_force[N] = 0.
                dihedral_tabulated = espressomd.interactions.TabulatedDihedral(
                    energy=tab_energy, force=tab_force)
                self.system.bonded_inter.add(dihedral_tabulated)
                self.system.part.clear()
                p0, p1, p2, p3 = self.system.part.add(pos=4 * [(0., 0., 0.)])
                p1.add_bond((dihedral_tabulated, p0, p2, p3))
                p1.pos = [5., 5., 5.]
                p2.pos = p1.pos + [1., 0., 0.]
                p0.pos = p1.pos + [0., 1., 0.]

                # use half the angular resolution to observe interpolation
                for i in range(2 * N - 1):
                    phi = i * d_phi / 2.
                    p3.pos = p2.pos + rotate_vector([0., 1., 0.], axis, phi)
                    self.system.integrator.run(recalc_forces=True, steps=0)

                    # Calculate expected forces and energies
                    j = i // 2
                    if i % 2 == 0:
                        E_ref = tab_energy[j]
                        _, forces_ref = dihedral_potential_and_forces(
                            dh_k, dh_n, dh_phi0, p0.pos, p1.pos, p2.pos, p3.pos)
                    else:
                        E_ref = (tab_energy[j] + tab_energy[j + 1]) / 2.0
                        forces_ref = None

                    self.check_values(E_ref, forces_ref)


if __name__ == '__main__':
    ut.main()
