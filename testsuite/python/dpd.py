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
import itertools

import espressomd
from espressomd.observables import DPDStress


@utx.skipIfMissingFeatures("DPD")
class DPDThermostat(ut.TestCase):

    """Test DPD dynamics."""

    s = espressomd.System(box_l=3 * [10.0])
    s.time_step = 0.01
    s.cell_system.skin = 0.4

    def setUp(self):
        np.random.seed(16)

    def tearDown(self):
        s = self.s
        s.part.clear()
        s.thermostat.turn_off()
        s.integrator.set_vv()

    def test_01__rng(self):
        """Test for RNG consistency."""
        def reset_particles():
            self.s.part.clear()
            p = self.s.part.add(pos=[0, 0, 0])
            _ = self.s.part.add(pos=[0, 0, 1])
            return p

        system = self.s

        kT = 2.3
        gamma = 1.5

        # No seed should throw exception
        with self.assertRaises(ValueError):
            system.thermostat.set_dpd(kT=kT)

        system.thermostat.set_dpd(kT=kT, seed=41)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)

        # run(0) does not increase the philox counter and should give the same
        # force
        p = reset_particles()
        system.integrator.run(0, recalc_forces=True)
        force0 = np.copy(p.f)
        system.integrator.run(0, recalc_forces=True)
        force1 = np.copy(p.f)
        np.testing.assert_almost_equal(force0, force1)

        # run(1) should give a different force
        p = reset_particles()
        system.integrator.run(1)
        force2 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force1, force2)))

        # Different seed should give a different force with same counter state
        # force2: dpd.rng_counter() = 1, dpd.rng_seed() = 41
        # force3: dpd.rng_counter() = 1, dpd.rng_seed() = 42
        p = reset_particles()
        system.integrator.run(0, recalc_forces=True)
        force2 = np.copy(p.f)
        system.thermostat.set_dpd(kT=kT, seed=42)
        system.integrator.run(0, recalc_forces=True)
        force3 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force2, force3)))

        # Same seed should not give the same force with different counter state
        p = reset_particles()
        system.thermostat.set_dpd(kT=kT, seed=42)
        system.integrator.run(1)
        force4 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force3, force4)))

        # Seed offset should not give the same force with a lag
        # force4: dpd.rng_counter() = 2, dpd.rng_seed() = 42
        # force5: dpd.rng_counter() = 3, dpd.rng_seed() = 41
        reset_particles()
        system.thermostat.set_dpd(kT=kT, seed=41)
        system.integrator.run(1)
        p = reset_particles()
        system.integrator.run(0, recalc_forces=True)
        force5 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force4, force5)))

    def test_const_weight_function(self):
        s = self.s
        kT = 0.
        gamma = 1.42
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.2,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.4)

        p0 = s.part.add(pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, .8, .3])
        p1 = s.part.add(pos=[3, 5, 5], type=0, v=v)

        s.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            np.testing.assert_array_equal(np.copy(f), [0., 0., 0.])

        # Only trans
        p1.pos = [5. - 1.3, 5, 5]

        s.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertLess(abs(p0.f[0]), 1e-16)
        np.testing.assert_allclose(
            np.copy(p0.f[1:2]), gamma * v[1:2], rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

        # Trans and parallel
        p1.pos = [5. - 1.1, 5, 5]

        s.integrator.run(0)

        np.testing.assert_allclose(
            np.copy(p0.f), gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

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

        p0 = s.part.add(pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, .8, .3])
        p1 = s.part.add(pos=[3, 5, 5], type=0, v=v)

        s.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            np.testing.assert_array_equal(np.copy(f), [0., 0., 0.])

        # Only trans
        p1.pos = [5. - 1.3, 5, 5]

        s.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertLess(abs(p0.f[0]), 1e-16)
        omega = calc_omega(1.3, 1.4)**2
        np.testing.assert_allclose(
            np.copy(p0.f[1:2]), omega * gamma * v[1:2], rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

        # Trans and parallel
        p1.pos = [5. - 1.1, 5, 5]

        s.integrator.run(0)

        omega = np.array([calc_omega(1.1, x)**2 for x in [1.2, 1.4, 1.4]])
        np.testing.assert_allclose(
            np.copy(p0.f), omega * gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

        # Trans and parallel 2nd point
        p1.pos = [5. - 0.5, 5, 5]

        s.integrator.run(0)

        omega = np.array([calc_omega(0.5, x)**2 for x in [1.2, 1.4, 1.4]])
        np.testing.assert_allclose(
            np.copy(p0.f), omega * gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

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

        p0 = s.part.add(pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, 0., 0.])
        p1 = s.part.add(pos=[3, 5, 5], type=0, v=v)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            np.testing.assert_array_equal(np.copy(f), [0., 0., 0.])

        # Place the particle at different positions to test the parabolic
        # weight function
        for dist in np.arange(0.1, 1.2, 50):

            p1.pos = [5. + dist, 5., 5.]
            s.integrator.run(0)
            omega = calc_omega(dist, r_cut)**2

            # The particle is moved along the x-direction. Hence, we are
            # testing the x element.
            np.testing.assert_allclose(
                np.copy(p0.f), omega * gamma * v, rtol=0, atol=1e-11)
            np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

    def test_ghosts_have_v(self):
        s = self.s

        r_cut = 1.5
        dx = 0.25 * r_cut
        ind_combinations = list(itertools.product([0, 1], [0, 1], [0, 1]))

        def f(i):
            if i == 0:
                return dx
            return 10. - dx

        # Put a particle in every corner
        for ind in ind_combinations:
            pos = [f(x) for x in ind]
            v = ind
            s.part.add(pos=pos, v=v)

        gamma = 1.0
        s.thermostat.set_dpd(kT=0.0, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=r_cut,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

        s.integrator.run(0)

        for p, ind in zip(s.part, ind_combinations):
            for i in ind:
                if ind[i] == 0:
                    sgn = 1
                else:
                    sgn = -1
                self.assertAlmostEqual(sgn * 4.0, p.f[i])

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

        np.testing.assert_array_almost_equal(np.copy(p.f), [1., 2., 3.])

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

        n_part = 200
        r_cut = 1.0
        gamma = 5.
        r_cut = 1.0

        s = self.s

        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=gamma / 2.0, trans_r_cut=r_cut)

        pos = s.box_l * np.random.random((n_part, 3))
        s.part.add(pos=pos)
        s.integrator.run(10)

        for kT in [0., 2.]:
            s.thermostat.set_dpd(kT=kT, seed=3)
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

            stress /= s.volume()

            dpd_stress = s.analysis.dpd_stress()

            dpd_obs = DPDStress()
            obs_stress = dpd_obs.calculate()

            np.testing.assert_array_almost_equal(np.copy(dpd_stress), stress)
            np.testing.assert_array_almost_equal(np.copy(obs_stress), stress)

    def test_momentum_conservation(self):
        r_cut = 1.0
        gamma = 5.
        r_cut = 2.9

        s = self.s
        s.thermostat.set_dpd(kT=1.3, seed=42)
        s.part.add(pos=((0, 0, 0), (0.1, 0.1, 0.1), (0.1, 0, 0)),
                   mass=(1, 2, 3))

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
