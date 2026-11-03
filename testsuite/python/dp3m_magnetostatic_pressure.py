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
import espressomd.magnetostatics


class pressureViaVolumeScaling:
    """
    Estimate a pressure contribution from isothermal volume scaling,
    without evaluating the virial. See "Efficient pressure estimation
    in molecular simulations without evaluating the virial" by
    Harismiadis, V. I., J. Vorholz, and A. Z. Panagiotopoulos (1996).
    Only valid for isotropic volume changes.

    ``energy_key`` selects which potential-energy contribution enters
    Delta U:

    * ``"total"`` reproduces the full mechanical pressure (matches
      ``pressure()['total']``) via the original Harismiadis identity,
      including the ideal-gas volume-rescaling term.
    * any other key (e.g. ``"dipolar"``) isolates a single additive
      energy contribution. Since dV/V is kept infinitesimal, the
      log-average linearizes to <Delta U>/kT, and <Delta U_total> is
      exactly the sum of the <Delta U_i> of independent contributions;
      this is the standard mean-force ("test-volume") estimator for a
      single contribution, so the ideal-gas term must NOT be added
      here (it would only apply to the total pressure).
    """

    def __init__(self, system, kT, energy_key="total"):
        self.system = system
        self.kT = kT
        self.energy_key = energy_key
        self.old_box_lengths = np.copy(system.box_l)
        self.old_volume = np.prod(self.old_box_lengths)
        dV_div_old_volume = 0.001
        self.dV = -dV_div_old_volume * self.old_volume
        self.new_volume = self.old_volume + self.dV
        self.new_box_l = (self.new_volume)**(1. / 3.)

        self.list_of_previous_values = []

    def _potential_energy(self):
        energy = self.system.analysis.energy()
        if self.energy_key == "total":
            return energy["total"] - energy["kinetic"]
        return energy[self.energy_key]

    def measure_pressure_via_volume_scaling(self):
        Epot_old = self._potential_energy()
        self.system.change_volume_and_rescale_particles(self.new_box_l, "xyz")
        self.system.integrator.run(0)
        Epot_new = self._potential_energy()
        self.system.change_volume_and_rescale_particles(
            self.old_box_lengths[0], "xyz")
        self.system.integrator.run(0)
        DeltaEpot = Epot_new - Epot_old

        if self.energy_key == "total":
            particle_number = len(self.system.part)
            current_value = (self.new_volume / self.old_volume)**particle_number * \
                np.exp(-DeltaEpot / self.kT)
        else:
            # mean-force estimator for a single additive energy
            # contribution: no ideal-gas Jacobian term, and no need
            # for the log/exp since dV is already infinitesimal
            current_value = -DeltaEpot / self.dV
        self.list_of_previous_values.append(current_value)

    def get_result(self):
        average_value = np.mean(self.list_of_previous_values)
        if self.energy_key == "total":
            return self.kT / self.dV * np.log(average_value)
        return average_value


@utx.skipIfMissingFeatures(["DP3M", "LENNARD_JONES"])
class VirialPressureConsistency(ut.TestCase):
    """
    Test the dipolar long-range (k-space) virial pressure against an
    analytical volume-scaling estimate (see
    :class:`pressureViaVolumeScaling`), and validate the full
    pressure tensor.

    An isotropic volume change alone can't validate the anisotropic
    pressure tensor: it's blind to off-diagonal terms, and passes
    even if the trace is split incorrectly among xx/yy/zz (e.g. the
    old ``diag(E, E, E) / 3`` placeholder). A strain-based check
    isn't available either, since DipolarP3M enforces a cubic box.
    :meth:`test_dp3m_pressure_tensor_rotation` covers this gap
    instead, using exact coordinate-relabeling identities that need
    no box deformation.
    """
    system = espressomd.System(box_l=[50, 50, 50])

    def setUp(self):
        np.random.seed(seed=1)
        self.system.time_step = 0.01
        self.kT = 0.5
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2**(1.0 / 6.0), shift="auto")
        num_part = 40
        mass = 1
        dip_magnitude = 1.3

        for _ in range(num_part):
            dip_dir = np.random.normal(size=3)
            dip_dir /= np.linalg.norm(dip_dir)
            self.system.part.add(
                pos=np.random.random(3) * self.system.box_l,
                dip=dip_magnitude * dip_dir,
                # keep dipole orientations fixed: the volume-scaling
                # estimator and the long-range pressure formula both
                # assume dipole moments are strain-invariant
                rotation=(False, False, False),
                v=np.sqrt(self.kT / mass) * np.random.normal(loc=[0, 0, 0]))

        #############################################################
        #  Warmup Integration                                       #
        #############################################################

        self.system.integrator.set_steepest_descent(
            f_max=0,
            gamma=0.001,
            max_displacement=0.01)

        # warmup
        energy = self.system.analysis.energy()["total"]
        print(f"minimization: {energy:.1f}")
        for _ in range(10):
            self.system.integrator.run(10)
            energy = self.system.analysis.energy()["total"]
            print(f"minimization: {energy:.1f}")
            if energy < 2 * num_part:
                break
        self.system.integrator.set_vv()
        self.system.thermostat.set_langevin(kT=self.kT, gamma=1.0, seed=42)
        # reset thermostat state and pin skin value to improve reproducibility
        self.system.thermostat.langevin.call_method(
            "override_philox_counter", counter=0)
        self.system.cell_system.skin = 1.4

    def tearDown(self):
        self.system.part.clear()
        self.system.magnetostatics.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def test_dp3m_pressure(self):
        self.system.magnetostatics.solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2., accuracy=1e-4, mesh=32, cao=6, r_cut=7.5, tune=True)
        skin = self.system.cell_system.tune_skin(
            min_skin=0.0, max_skin=2.5, tol=0.05, int_steps=100)
        print(f"Tuned skin: {skin}")

        num_samples = 25
        pressures_via_virial = {"total": [], "dipolar": []}
        pressure_via_volume_scaling_total = pressureViaVolumeScaling(
            self.system, self.kT, energy_key="total")
        pressure_via_volume_scaling_dipolar = pressureViaVolumeScaling(
            self.system, self.kT, energy_key="dipolar")
        for _ in range(num_samples):
            self.system.integrator.run(50)
            pressure = self.system.analysis.pressure()
            pressures_via_virial["total"].append(pressure["total"])
            pressures_via_virial["dipolar"].append(pressure["dipolar"])
            pressure_via_volume_scaling_total.measure_pressure_via_volume_scaling()
            pressure_via_volume_scaling_dipolar.measure_pressure_via_volume_scaling()

        # deviation should be below 5%
        pressure_virial_total = np.mean(pressures_via_virial["total"])
        abs_deviation_total = 100 * abs(
            pressure_virial_total /
            pressure_via_volume_scaling_total.get_result() - 1.0)
        np.testing.assert_array_less(abs_deviation_total, 5.0)

        # isolated check of the dipolar long-range pressure term alone.
        pressure_virial_dipolar = np.mean(pressures_via_virial["dipolar"])
        pressure_scaling_dipolar = pressure_via_volume_scaling_dipolar.get_result()
        np.testing.assert_allclose(
            pressure_virial_dipolar, pressure_scaling_dipolar,
            rtol=0.1, atol=1e-6)

    def test_dp3m_pressure_tensor_rotation(self):
        """
        Check that the k-space dipolar pressure tensor transforms as a
        genuine rank-2 tensor under two exact coordinate relabelings: a
        90 degree rotation about the z-axis, and a cyclic permutation of
        the axes. Both only relabel coordinates (no box deformation),
        which matters because DipolarP3M enforces a cubic box.

        The k-space tensor is generally asymmetric, so the expected
        rotated/permuted tensors are built with the general covariant
        transformation law ``sigma' = R sigma R^T`` which is valid for
        any rank-2 tensor.
        """
        self.system.magnetostatics.solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2., accuracy=1e-4, mesh=32, cao=6, r_cut=7.5, tune=True)
        self.system.integrator.run(0)

        box_l = self.system.box_l[0]
        pos = self.system.part.all().pos
        dip = self.system.part.all().dip

        def get_kspace_tensor():
            return self.system.analysis.pressure_tensor()[("dipolar", 1)]

        def set_configuration(new_pos, new_dip):
            self.system.part.all().pos = new_pos
            self.system.part.all().dip = new_dip
            self.system.integrator.run(0)

        sigma = get_kspace_tensor()

        # sanity check: the diagonal entries must not all be equal, and
        # the off-diagonal entries must not all vanish -- otherwise the
        # diagonal or off-diagonal half of the symmetry checks below
        # would hold vacuously (as they would for the old isotropic
        # diag(E, E, E) / 3 placeholder).
        scale = np.max(np.abs(sigma))
        diagonal_spread = np.std(np.diag(sigma))
        offdiagonal_scale = np.max(
            np.abs([sigma[0, 1], sigma[0, 2], sigma[1, 2]]))
        antisymmetric_scale = np.max(np.abs(sigma - sigma.T))
        self.assertGreater(diagonal_spread, 1e-3 * scale)
        self.assertGreater(offdiagonal_scale, 1e-3 * scale)
        self.assertGreater(antisymmetric_scale, 1e-3 * scale)

        # 90 degree rotation about the z-axis through the box center
        rot_pos = np.copy(pos)
        rot_pos[:, 0] = np.mod(box_l - pos[:, 1], box_l)
        rot_pos[:, 1] = np.mod(pos[:, 0], box_l)
        rot_dip = np.copy(dip)
        rot_dip[:, 0] = -dip[:, 1]
        rot_dip[:, 1] = dip[:, 0]
        set_configuration(rot_pos, rot_dip)
        sigma_rotated = get_kspace_tensor()

        R_rot = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        expected_rotated = R_rot @ sigma @ R_rot.T
        np.testing.assert_allclose(
            sigma_rotated, expected_rotated, atol=1e-9 * scale)

        # cyclic permutation of the axes (about the original, unrotated
        # configuration, restored here for an independent check)
        set_configuration(pos, dip)
        perm_pos = pos[:, [1, 2, 0]]
        perm_dip = dip[:, [1, 2, 0]]
        set_configuration(perm_pos, perm_dip)
        sigma_permuted = get_kspace_tensor()

        R_perm = np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])
        expected_permuted = R_perm @ sigma @ R_perm.T
        np.testing.assert_allclose(
            sigma_permuted, expected_permuted, atol=1e-9 * scale)

    def test_dp3m_pressure_tensor_vs_continuum_ewald(self):
        """
        Validate the dipolar reciprocal-space pressure tensor against
        an independent, from-scratch calculation: a direct sum over
        explicit wavevectors in Python, with no FFT, no mesh, and no
        shared code with P3M.

        This complements :func:`test_dp3m_pressure_tensor_rotation`,
        which only checks that the tensor transforms correctly under
        rotation -- a check that a uniformly wrong prefactor (e.g. a
        missing factor of 2) would still pass. Comparing against an
        independently derived formula catches that class of error.

        The symmetric part of the reference matches eq. (46) of
        :cite:`aguado03a` (up to a factor of 2, from summing here over
        all ``k != 0`` instead of their half-space ``h>0``). That
        paper only reports the symmetric form; the antisymmetric part
        is derived here independently by differentiating the
        reciprocal energy under a general strain (matching the
        derivation in ``long_range_pressure()`` from DP3M).
        """
        # a tighter accuracy than the other tests in this module (which
        # use 1e-4) is requested here, letting mesh/cao auto-tune to a
        # finer mesh: this test compares against the *unmeshed*
        # continuum sum, so P3M's own mesh-discretization error must be
        # small enough not to dominate the comparison tolerance below
        self.system.magnetostatics.solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2., accuracy=1e-6, r_cut=7.5, tune=True)
        self.system.integrator.run(0)

        solver = self.system.magnetostatics.solver
        alpha = solver.alpha
        prefactor = solver.prefactor
        box_l = self.system.box_l[0]
        volume = box_l**3
        pos = self.system.part.all().pos
        dip = self.system.part.all().dip

        sigma_p3m = self.system.analysis.pressure_tensor()[("dipolar", 1)]

        # brute-force continuum sum over wavevectors k = (2 pi / L) n;
        # n_max is chosen generously large so that exp(-k^2/(4 alpha^2))
        # has decayed to a negligible size well within the cutoff
        n_max = 40
        ns = np.arange(-n_max, n_max + 1)
        nx, ny, nz = np.meshgrid(ns, ns, ns, indexing="ij")
        nvecs = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=-1)
        nvecs = nvecs[np.any(nvecs != 0, axis=1)]
        kvecs = (2. * np.pi / box_l) * nvecs
        k2 = np.sum(kvecs**2, axis=1)
        converged = k2 / (4. * alpha**2) < 40.
        kvecs = kvecs[converged]
        k2 = k2[converged]

        # M_a(k) = sum_j mu_j,a exp(i k.r_j); Q(k) = M(k).k
        phases = np.exp(1j * (kvecs @ pos.T))
        M = phases @ dip
        Q = np.einsum("ka,ka->k", kvecs, M)

        g = np.exp(-k2 / (4. * alpha**2)) / k2
        vterm = -2. * (1. / k2 + 1. / (4. * alpha**2))
        cell_energy = np.abs(Q)**2

        # full (generally asymmetric) tensor: Pi_ab = g * [(delta_ab +
        # vterm * k_a * k_b) * cell_energy + 2 * k_b * Re(M_a(k) Q(k)^*)];
        # its symmetrized form (average with the transpose) is eq. (46)
        # in Aguado & Madden. Note the cross term's indices are swapped
        # relative to the strain probe that produces it (probing
        # epsilon_ab yields a term proportional to k_a Re[M_b Q*], which
        # belongs at tensor position (b, a) -- see the class-level
        # comment on ``long_range_pressure()`` in ``dp3m_heffte.impl.hpp`` for
        # the derivation): a pair virial is conventionally r_a F_b, but
        # differentiating the reciprocal energy under the strain probe
        # H(eps) = L*I + eps*E_ab yields r_b F_a.
        reference = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                diag = cell_energy if a == b else 0.
                envelope = cell_energy * vterm * kvecs[:, a] * kvecs[:, b]
                re_Ma_Qc = (M[:, a] * np.conj(Q)).real
                cross = 2. * kvecs[:, b] * re_Ma_Qc
                reference[a, b] = np.sum(g * (diag + envelope + cross))
        reference *= (2. * np.pi / volume**2) * prefactor

        scale = np.max(np.abs(reference))
        np.testing.assert_allclose(sigma_p3m, reference, atol=0.01 * scale)

    def test_dp3m_pressure_tensor_vs_direct_sum(self):
        """
        Check the total dipolar pressure tensor against an
        independent solver,
        :class:`~espressomd.magnetostatics.DipolarDirectSum`,
        exercising its C++ periodic-image summation end to end --
        otherwise untested for the pressure tensor.

        DP3M and DipolarDirectSum resolve the same conditionally
        convergent periodic lattice sum in different ways: DP3M
        assumes conducting ("tin-foil") boundary conditions, while
        DirectSum sums a growing *cube* of periodic images with no
        such correction. The two conventions differ by a term
        proportional to the square of the net dipole moment (de
        Leeuw, Perram & Smith, 1980), so this test uses opposite-dipole
        pairs to make that moment exactly zero and remove the
        ambiguity.

        This also compares the *antisymmetric* (torque-related) part of
        the pressure tensor. Unlike the symmetric part, it has no
        boundary/shape ambiguity to begin with: the de Leeuw-Perram-Smith
        term is an isotropic function of the net dipole moment alone (no
        dependence on individual particle positions), so it can only ever
        contribute to the trace, never to any off-diagonal component --
        symmetric or antisymmetric. So DP3M and DirectSum must agree on
        the antisymmetric part directly, independent of net dipole moment.
        """
        system = self.system
        system.part.clear()
        system.thermostat.turn_off()

        np.random.seed(seed=2)
        n_pairs = 8
        box_l = system.box_l[0]
        pos = np.random.random((2 * n_pairs, 3)) * box_l
        dip = np.random.normal(size=(n_pairs, 3))
        dip = np.concatenate([dip, -dip])
        system.part.add(pos=pos, dip=dip,
                        rotation=(2 * n_pairs) * [(False, False, False)])

        # by construction: pairs of opposite dipole moments
        net_dipole_moment = np.sum(system.part.all().dip, axis=0)
        np.testing.assert_allclose(net_dipole_moment, 0., atol=1e-9)

        system.magnetostatics.solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2., accuracy=1e-6, r_cut=7.5, tune=True)
        system.integrator.run(0)
        pressure_p3m = system.analysis.pressure_tensor()["dipolar"]

        system.magnetostatics.solver = espressomd.magnetostatics.DipolarDirectSum(
            prefactor=2., n_replicas=6)
        system.integrator.run(0)
        pressure_dds = system.analysis.pressure_tensor()["dipolar"]

        sym_p3m = (pressure_p3m + pressure_p3m.T) / 2.
        sym_dds = (pressure_dds + pressure_dds.T) / 2.
        antisym_p3m = (pressure_p3m - pressure_p3m.T) / 2.
        antisym_dds = (pressure_dds - pressure_dds.T) / 2.
        scale = np.max(np.abs(sym_p3m))

        # sanity check: the antisymmetric part must not vanish, otherwise
        # the comparison below would hold vacuously
        self.assertGreater(np.max(np.abs(antisym_p3m)), 1e-3 * scale)

        np.testing.assert_allclose(
            np.trace(pressure_dds), np.trace(pressure_p3m), rtol=1e-2)
        np.testing.assert_allclose(sym_dds, sym_p3m, atol=0.02 * scale)
        np.testing.assert_allclose(antisym_dds, antisym_p3m, atol=0.02 * scale)


if __name__ == "__main__":
    ut.main()
