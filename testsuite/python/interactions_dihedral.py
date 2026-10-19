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


def rotate_vector(v, k, phi):
    """Rotates vector v around unit vector k by angle phi.
    Uses Rodrigues' rotation formula."""
    vrot = np.array(v) * np.cos(phi) + np.cross(k, v) * \
        np.sin(phi) + np.array(k) * np.dot(k, v) * (1.0 - np.cos(phi))
    return vrot


def dihedral_angle(p1, p2, p3, p4):
    """Dihedral angle of a particle quadruple, in the range [0, 2 pi).
    """
    v12, v23, v34 = p2 - p1, p3 - p2, p4 - p3
    n1 = np.cross(v12, v23)
    n2 = np.cross(v23, v34)
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    cosphi = np.clip(n1 @ n2, -1., 1.)
    # acos cannot distinguish phi from 2pi-phi
    phi = np.arccos(cosphi)
    # the sign of n1.v34 helps find a unique phi from 0 to 2pi (full range)
    if (n1 @ v34) < 0.:
        phi = 2. * np.pi - phi
    return phi


def dihedral_angle_gradients(p1, p2, p3, p4):
    """Gradients of the dihedral angle w.r.t. the four positions.

    Taken from Blondel and Karplus, J. Comput. Chem. 17, 1132 (1996), eq. 27.
    Unlike expressions derived by differentiating cos(phi), these carry no
    1/sin(phi) factor and are therefore valid at phi = 0 and phi = pi too.
    """
    v12, v23, v34 = p2 - p1, p3 - p2, p4 - p3
    # Blondel's F = -v12, G = -v23, H = v34, so A = F x G and B = H x G are
    A = np.cross(v12, v23)
    B = np.cross(v23, v34)
    A_sqr, B_sqr = np.dot(A, A), np.dot(B, B)
    l_G = np.linalg.norm(v23)

    cA, cB = l_G / A_sqr, l_G / B_sqr
    dA = np.dot(v12, v23) / (A_sqr * l_G)    # (F.G) / (A^2 |G|)
    dB = -np.dot(v23, v34) / (B_sqr * l_G)   # (H.G) / (B^2 |G|)

    return (-cA * A,
            (cA + dA) * A - dB * B,
            (dB - cB) * B - dA * A,
            cB * B)


def dihedral_potential_and_forces(k, n, phase, p1, p2, p3, p4):
    """
    Reference potential and forces for a dihedral angle.

    The force follows from the plain chain rule
    ``f_i = -(dV/dphi) * dphi/dr_i``, with the gradients from
    :func:`dihedral_angle_gradients`. This reference is valid for every angle,
    including phi = 0 and phi = pi.
    """
    phi = dihedral_angle(p1, p2, p3, p4)
    dV_dphi = k * n * np.sin(n * phi - phase)
    forces = tuple(-dV_dphi * grad
                   for grad in dihedral_angle_gradients(p1, p2, p3, p4))
    potential = k * (1 - np.cos(n * phi - phase))
    return (potential, forces)


def forces_from_energy_gradient(energy_of_phi, p1, p2, p3, p4, h=2e-4):
    """Forces obtained by numerically differentiating the energy.

    It is independent of any closed-form force expression: it only ever
    evaluates the energy as a function of the dihedral angle. A force formula
    that is wrong in both the core and the analytic reference 
    (def dihedral_potential_and_forces()) above would still
    be caught here.

    The step size h=2e-4 is chosen such that it is not too small for acos(cos phi)
    operation to loose precision or nor too large such that finite difference truncation
    error is accumulated.
    """
    positions = [np.copy(p) for p in (p1, p2, p3, p4)]
    forces = []
    for atom in range(4):
        force = np.zeros(3)
        for axis in range(3):
            shifted = [np.copy(p) for p in positions]
            shifted[atom][axis] += h
            e_plus = energy_of_phi(dihedral_angle(*shifted))
            shifted[atom][axis] -= 2. * h
            e_minus = energy_of_phi(dihedral_angle(*shifted))
            force[axis] = -(e_plus - e_minus) / (2. * h)
        forces.append(force)
    return tuple(forces)


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

    def check_collinear_particles(self):
        """Three collinear particles leave the dihedral angle undefined: no
        plane normal can be constructed to measure the angle from, and (as
        derived by Blondel and Karplus (1996)) the force gradient genuinely diverges as any
        one of the three approaches the line through the other two.

        Matching LAMMPS's convention for this degenerate case, the core
        raises a runtime warning and reports a zero energy/force/pressure
        for this bond, rather than treating it as a broken bond.
        """
        p0 = self.system.part.by_id(0)
        p1 = self.system.part.by_id(1)
        p2 = self.system.part.by_id(2)
        p3 = self.system.part.by_id(3)
        p1.pos = [5., 5., 5.]
        p2.pos = p1.pos + [1., 0., 0.]
        p3.pos = p1.pos + [2., 0., 0.]
        p0.pos = p1.pos + [-1., 0., 0.]
        self.assertEqual(self.system.analysis.energy()["bonded"], 0.)
        self.system.integrator.run(steps=0, recalc_forces=True)
        f0, f1, f2, f3 = self.system.part.all().f
        np.testing.assert_array_equal(np.copy(f0), np.zeros(3))
        np.testing.assert_array_equal(np.copy(f1), np.zeros(3))
        np.testing.assert_array_equal(np.copy(f2), np.zeros(3))
        np.testing.assert_array_equal(np.copy(f3), np.zeros(3))

    def check_undefined_angle(self):
        """A zero-length bond vector (coincident particles) is degenerate
        for the same underlying reason as three collinear particles: no
        plane normal can be constructed. Handled identically: a runtime
        warning, zero energy/force, no exception.
        """
        p0 = self.system.part.by_id(0)
        p1 = self.system.part.by_id(1)
        p0.pos = p1.pos
        self.assertEqual(self.system.analysis.energy()["bonded"], 0.)
        self.system.integrator.run(steps=0, recalc_forces=True)
        f0, f1, f2, f3 = self.system.part.all().f
        np.testing.assert_array_equal(np.copy(f0), np.zeros(3))
        np.testing.assert_array_equal(np.copy(f1), np.zeros(3))
        np.testing.assert_array_equal(np.copy(f2), np.zeros(3))
        np.testing.assert_array_equal(np.copy(f3), np.zeros(3))

    def check_pressure_tensor(self, tol=1e-12):
        p0, p1, p2, p3 = self.system.part.all()
        # p1 is the bond owner (reference particle)
        # P_ij = 1/V * sum F_{1,k}_i * r_{1,k}_j
        p_tensor_ref = (
            np.outer(np.copy(p0.f), self.system.distance_vec(p1, p0))
            + np.outer(np.copy(p2.f), self.system.distance_vec(p1, p2))
            + np.outer(np.copy(p3.f), self.system.distance_vec(p1, p3))
        ) / self.system.volume()
        p_tensor_sim = self.system.analysis.pressure_tensor()["bonded"]
        np.testing.assert_allclose(p_tensor_sim,
                                   p_tensor_ref,
                                   atol=tol)
        np.testing.assert_allclose(self.system.analysis.pressure()["bonded"],
                                   0.0,
                                   atol=tol)
        # consistency: trace / 3 == scalar pressure
        np.testing.assert_allclose(np.trace(p_tensor_sim) / 3.,
                                   self.system.analysis.pressure()["bonded"],
                                   atol=tol)
        # symmetry: dihedral angle is rotationally invariant, so
        # sum_i r_i x F_i = 0, which implies p_ab = p_ba
        np.testing.assert_allclose(p_tensor_sim,
                                   p_tensor_sim.T,
                                   atol=tol)

    def check_no_singularity_at_0_and_pi(self, k, n, phase):
        """The force at phi = 0 and phi = pi must be finite and correct.

        When the phase is not a multiple of pi, dV/dphi does not vanish at
        these two angles, so the force there is non-zero. Formulations that
        differentiate cos(phi) pick up a 1/sin(phi) factor which is unbounded
        there; this asserts that no such artefact is left.
        """
        axis = np.array([1., 0., 0.])
        p0, p1, p2, p3 = self.system.part.all()
        p1.pos = [5., 5., 5.]
        p2.pos = p1.pos + [1., 0., 0.]
        p0.pos = p1.pos + [0., 1., 0.]
        for phi in [0., np.pi]:
            with self.subTest(phi=phi):
                p3.pos = p2.pos + rotate_vector([0., 1., 0.], axis, phi)
                self.system.integrator.run(recalc_forces=True, steps=0)
                # p0 and p3 both sit one length unit off the p1-p2 axis, so
                # |dphi/dr| = 1 for them and hence |f| = |dV/dphi|
                f_ref = abs(k * n * np.sin(n * phi - phase))
                for p in (p0, p3):
                    np.testing.assert_allclose(
                        np.linalg.norm(np.copy(p.f)), f_ref, rtol=1e-10)
                E_ref, forces_ref = dihedral_potential_and_forces(
                    k, n, phase, p0.pos, p1.pos, p2.pos, p3.pos)
                self.check_values(E_ref, forces_ref)

    def test_forces_match_energy_gradient(self):
        """Cross-check the core forces against a numerical energy gradient.

        Independent of any closed-form force expression, so an error shared by
        the core implementation of Blondel and Karplus (1996)
        in this file would still be caught.
        """
        axis = np.array([1., 0., 0.])
        k, n, phase = 2., 2, np.pi / 3
        dihedral = espressomd.interactions.Dihedral(
            bend=k, mult=n, phase=phase)
        self.system.bonded_inter.add(dihedral)
        self.system.part.clear()
        p0, p1, p2, p3 = self.system.part.add(pos=4 * [(0., 0., 0.)])
        p1.add_bond((dihedral, p0, p2, p3))
        p1.pos = [5., 5., 5.]
        p2.pos = p1.pos + [1., 0., 0.]
        p0.pos = p1.pos + [0., 1., 0.]

        def energy_of_phi(phi):
            return k * (1. - np.cos(n * phi - phase))

        for phi in [0., 1e-6, 0.4, np.pi / 2., np.pi - 1e-6, np.pi, 4.,
                    2. * np.pi - 0.3]:
            with self.subTest(phi=phi):
                p3.pos = p2.pos + rotate_vector([0., 1., 0.], axis, phi)
                self.system.integrator.run(recalc_forces=True, steps=0)
                forces_ref = forces_from_energy_gradient(
                    energy_of_phi, np.copy(p0.pos), np.copy(p1.pos),
                    np.copy(p2.pos), np.copy(p3.pos))
                for p, f_ref in zip((p0, p1, p2, p3), forces_ref):
                    np.testing.assert_allclose(np.copy(p.f), f_ref, atol=1e-6)

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
                    self.check_pressure_tensor()

                self.check_no_singularity_at_0_and_pi(dh_k, dh_n, dh_phi0)

        self.check_undefined_angle()
        self.check_collinear_particles()

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
                # tabulated values for the range [0, 2*pi]; the force
                # table holds -dV/dphi, which is bounded everywhere
                phi = d_phi * np.arange(N + 1)
                tab_energy = dh_k * (1. - np.cos(dh_n * phi - dh_phi0))
                tab_force = -dh_k * dh_n * np.sin(dh_n * phi - dh_phi0)
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
                    self.check_pressure_tensor()

                self.check_no_singularity_at_0_and_pi(dh_k, dh_n, dh_phi0)

        self.check_undefined_angle()
        self.check_collinear_particles()


if __name__ == '__main__':
    ut.main()
