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
import numpy as np
import unittest as ut
import unittest_decorators as utx
from itertools import product

import espressomd
from espressomd.observables import DPDStress
from tests_common import single_component_maxwell


@utx.skipIfMissingFeatures("DPD")
class DPDThermostat(ut.TestCase):

    """Tests the velocity distribution created by the dpd thermostat against
       the single component Maxwell distribution."""

    s = espressomd.System(box_l=3 * [10.0])
    s.time_step = 0.01
    s.cell_system.skin = 0.4

    def setUp(self):
        np.random.seed(16)

    def tearDown(self):
        s = self.s
        s.part.clear()
        s.thermostat.turn_off()

    def check_velocity_distribution(self, vel, minmax, n_bins, error_tol, kT):
        """check the recorded particle distributions in velocity against a
           histogram with n_bins bins. Drop velocities outside minmax. Check
           individual histogram bins up to an accuracy of error_tol against
           the analytical result for kT."""
        for i in range(3):
            hist = np.histogram(
                vel[:, i], range=(-minmax, minmax), bins=n_bins, density=False)
            data = hist[0] / float(vel.shape[0])
            bins = hist[1]
            for j in range(n_bins):
                found = data[j]
                expected = single_component_maxwell(bins[j], bins[j + 1], kT)
                self.assertAlmostEqual(found, expected, delta=error_tol)

    def test_aa_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertAlmostEqual(
            single_component_maxwell(-10, 10, 4.), 1., delta=1E-4)

    def check_total_zero(self):
        v_total = np.sum(self.s.part[:].v, axis=0)
        np.testing.assert_allclose(v_total, np.zeros(3), atol=1e-11)

    def single(self, with_langevin=False):
        """Test velocity distribution of a dpd fluid with a single type."""
        N = 500
        s = self.s
        s.part.add(pos=s.box_l * np.random.random((N, 3)))
        kT = 2.3
        gamma = 1.5
        if with_langevin:
            s.thermostat.set_langevin(kT=kT, gamma=gamma, seed=41)
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(100)
        loops = 250
        v_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i * N:(i + 1) * N, :] = s.part[:].v
        v_minmax = 5
        bins = 5
        error_tol = 0.01
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)

        if not with_langevin:
            self.check_total_zero()

    def test_single(self):
        self.single()

    def test_single_with_langevin(self):
        self.single(True)

    def test_binary(self):
        """Test velocity distribution of binary dpd fluid"""
        N = 200
        s = self.s
        s.part.add(pos=s.box_l * np.random.random((N // 2, 3)),
                   type=N // 2 * [0])
        s.part.add(pos=s.box_l * np.random.random((N // 2, 3)),
                   type=N // 2 * [1])
        kT = 2.3
        gamma = 1.5
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.0,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.0)
        s.non_bonded_inter[1, 1].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.0,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.0)
        s.non_bonded_inter[0, 1].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(100)
        loops = 400
        v_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i * N:(i + 1) * N, :] = s.part[:].v
        v_minmax = 5
        bins = 5
        error_tol = 0.01
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)
        self.check_total_zero()

    def test_disable(self):
        N = 200
        s = self.s
        s.time_step = 0.01
        s.part.add(pos=s.box_l * np.random.random((N, 3)))
        kT = 2.3
        gamma = 1.5
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)

        s.integrator.run(10)

        s.thermostat.turn_off()

        # Reset velocities
        s.part[:].v = [1., 2., 3.]

        s.integrator.run(10)

        # Check that there was neither noise nor friction
        for v in s.part[:].v:
            for i in range(3):
                self.assertEqual(v[i], float(i + 1))

        # Turn back on
        s.thermostat.set_dpd(kT=kT, seed=42)

        # Reset velocities for faster convergence
        s.part[:].v = [0., 0., 0.]

        # Equilibrate
        s.integrator.run(250)

        loops = 250
        v_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i * N:(i + 1) * N, :] = s.part[:].v
        v_minmax = 5
        bins = 5
        error_tol = 0.012
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)

    def test_const_weight_function(self):
        s = self.s
        kT = 0.
        gamma = 1.42
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.2,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.4)

        s.part.add(id=0, pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, .8, .3])
        s.part.add(id=1, pos=[3, 5, 5], type=0, v=v)

        s.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            np.testing.assert_array_equal(f, [0., 0., 0.])

        # Only trans
        s.part[1].pos = [5. - 1.3, 5, 5]

        s.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertLess(abs(s.part[0].f[0]), 1e-16)
        np.testing.assert_allclose(
            np.copy(s.part[0].f[1:2]), gamma * v[1:2], rtol=0, atol=1e-11)
        np.testing.assert_array_equal(
            np.copy(s.part[0].f), -np.copy(s.part[1].f))

        # Trans and parallel
        s.part[1].pos = [5. - 1.1, 5, 5]

        s.integrator.run(0)

        np.testing.assert_allclose(
            np.copy(s.part[0].f), gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(
            np.copy(s.part[0].f), -np.copy(s.part[1].f))

    def test_linear_weight_function(self):
        s = self.s
        kT = 0.
        gamma = 1.42
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=1.2,
            trans_weight_function=1, trans_gamma=gamma, trans_r_cut=1.4)

        def calc_omega(dist, r_cut):
            return 1. - dist / r_cut

        s.part.add(id=0, pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, .8, .3])
        s.part.add(id=1, pos=[3, 5, 5], type=0, v=v)

        s.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            np.testing.assert_array_equal(f, [0., 0., 0.])

        # Only trans
        s.part[1].pos = [5. - 1.3, 5, 5]

        s.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertLess(abs(s.part[0].f[0]), 1e-16)
        omega = calc_omega(1.3, 1.4)**2
        np.testing.assert_allclose(
            np.copy(s.part[0].f[1:2]), omega * gamma * v[1:2], rtol=0, atol=1e-11)
        np.testing.assert_array_equal(
            np.copy(s.part[0].f), -np.copy(s.part[1].f))

        # Trans and parallel
        s.part[1].pos = [5. - 1.1, 5, 5]

        s.integrator.run(0)

        omega = np.array([calc_omega(1.1, x)**2 for x in [1.2, 1.4, 1.4]])
        np.testing.assert_allclose(
            np.copy(s.part[0].f), omega * gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(
            np.copy(s.part[0].f), -np.copy(s.part[1].f))

        # Trans and parallel 2nd point
        s.part[1].pos = [5. - 0.5, 5, 5]

        s.integrator.run(0)

        omega = np.array([calc_omega(0.5, x)**2 for x in [1.2, 1.4, 1.4]])
        np.testing.assert_allclose(
            np.copy(s.part[0].f), omega * gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(
            np.copy(s.part[0].f), -np.copy(s.part[1].f))

    def test_parabolic_weight_function(self):
        s = self.s
        kT = 0.
        gamma = 1.42
        kappa = 7.8
        r_cut = 1.2
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, k=kappa, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=0.0, trans_r_cut=0.0)

        def calc_omega(dist, r_cut):
            return (1. - (dist / r_cut) ** kappa)

        s.part.add(id=0, pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, 0., 0.])
        s.part.add(id=1, pos=[3, 5, 5], type=0, v=v)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            np.testing.assert_array_equal(f, [0., 0., 0.])

        # Place the particle at different positions to test the parabolic
        # weight function
        for dist in np.arange(0.1, 1.2, 50):

            s.part[1].pos = [5. + dist, 5., 5.]
            s.integrator.run(0)
            omega = calc_omega(dist, r_cut)**2

            # The particle is moved along the x-direction. Hence, we are
            # testing the x element.
            np.testing.assert_allclose(
                np.copy(s.part[0].f), omega * gamma * v, rtol=0, atol=1e-11)
            np.testing.assert_array_equal(
                np.copy(s.part[0].f), -np.copy(s.part[1].f))

    def test_ghosts_have_v(self):
        s = self.s

        r_cut = 1.5
        dx = 0.25 * r_cut

        def f(i):
            if i == 0:
                return dx
            return 10. - dx

        # Put a particle in every corner
        for ind in product([0, 1], [0, 1], [0, 1]):
            pos = [f(x) for x in ind]
            v = ind
            s.part.add(pos=pos, v=v)

        gamma = 1.0
        s.thermostat.set_dpd(kT=0.0, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=r_cut,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

        s.integrator.run(0)

        id = 0
        for ind in product([0, 1], [0, 1], [0, 1]):
            for i in ind:
                if ind[i] == 0:
                    sgn = 1
                else:
                    sgn = -1
                self.assertAlmostEqual(sgn * 4.0, s.part[id].f[i])
            id += 1

    def test_constraint(self):
        import espressomd.shapes

        s = self.s

        s.constraints.add(shape=espressomd.shapes.Wall(
            dist=0, normal=[1, 0, 0]), particle_type=0, particle_velocity=[1, 2, 3])

        s.thermostat.set_dpd(kT=0.0, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=1., r_cut=1.0,
            trans_weight_function=0, trans_gamma=1., trans_r_cut=1.0)

        p = s.part.add(pos=[0.5, 0, 0], type=0, v=[0, 0, 0])

        s.integrator.run(0)

        self.assertAlmostEqual(p.f[0], 1.)
        self.assertAlmostEqual(p.f[1], 2.)
        self.assertAlmostEqual(p.f[2], 3.)

        for c in s.constraints:
            s.constraints.remove(c)

    def test_dpd_stress(self):

        def calc_omega(dist):
            return (1. / dist - 1. / r_cut) ** 2.0

        def diss_force_1(dist, vel_diff):
            f = np.zeros(3)
            vel12dotd12 = 0.
            dist_norm = np.linalg.norm(dist)

            for d in range(3):
                vel12dotd12 += vel_diff[d] * dist[d]

            friction = gamma * calc_omega(dist_norm) * vel12dotd12

            for d in range(3):
                f[d] -= (dist[d] * friction)

            return f

        def diss_force_2(dist, vel_diff):
            dist_norm = np.linalg.norm(dist)
            mat = np.identity(3) * (dist_norm**2.0)

            f = np.zeros(3)

            for d1 in range(3):
                for d2 in range(3):
                    mat[d1, d2] -= dist[d1] * dist[d2]

            for d1 in range(3):
                for d2 in range(3):
                    f[d1] += mat[d1, d2] * vel_diff[d2]
                f[d1] *= - 1.0 * gamma / 2.0 * calc_omega(dist_norm)

            return f

        def calc_stress(dist, vel_diff):
            force_pair = diss_force_1(dist, vel_diff) +\
                diss_force_2(dist, vel_diff)
            stress_pair = np.outer(dist, force_pair)
            return stress_pair

        n_part = 1000
        r_cut = 1.0
        gamma = 5.
        r_cut = 1.0

        s = self.s
        s.part.clear()

        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=gamma / 2.0, trans_r_cut=r_cut)

        pos = s.box_l * np.random.random((n_part, 3))
        s.part.add(pos=pos)
        s.integrator.run(10)

        for kT in [0., 2.]:
            s.thermostat.set_dpd(kT=kT)
            # run 1 integration step to get velocities
            s.part[:].v = np.zeros((n_part, 3))
            s.integrator.run(steps=1)

            pairs = s.part.pairs()

            stress = np.zeros([3, 3])

            for pair in pairs:
                dist = s.distance_vec(pair[0], pair[1])
                if np.linalg.norm(dist) < r_cut:
                    vel_diff = pair[1].v - pair[0].v
                    stress += calc_stress(dist, vel_diff)

            stress /= np.prod(np.copy(s.box_l))

            dpd_stress = s.analysis.dpd_stress()

            dpd_obs = DPDStress()
            obs_stress = np.array(dpd_obs.calculate()).reshape((3, 3))

            np.testing.assert_array_almost_equal(np.copy(dpd_stress), stress)
            np.testing.assert_array_almost_equal(np.copy(obs_stress), stress)

    def test_momentum_conservation(self):
        r_cut = 1.0
        gamma = 5.
        r_cut = 2.9

        s = self.s
        s.thermostat.set_dpd(kT=1.3, seed=42)
        s.part.clear()
        s.part.add(pos=((0, 0, 0), (0.1, 0.1, 0.1),
                        (0.1, 0, 0)), mass=(1, 2, 3))

        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=gamma / 2.0, trans_r_cut=r_cut)
        momentum = np.matmul(s.part[:].v.T, s.part[:].mass)
        for _ in range(10):
            s.integrator.run(25)
            np.testing.assert_almost_equal(np.zeros((3,)), np.sum(s.part[:].f))
            np.testing.assert_allclose(
                np.matmul(s.part[:].v.T, s.part[:].mass), momentum, atol=1E-12)


if __name__ == "__main__":
    ut.main()
