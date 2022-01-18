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
import espressomd.observables
import espressomd.constraints
import espressomd.shapes


@utx.skipIfMissingFeatures("DPD")
class DPDThermostat(ut.TestCase):

    """Test DPD dynamics."""

    system = espressomd.System(box_l=3 * [10.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        np.random.seed(16)

    def tearDown(self):
        self.system.part.clear()
        self.system.constraints.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def test_01__rng(self):
        """Test for RNG consistency."""
        def reset_particles():
            self.system.part.clear()
            p = self.system.part.add(pos=[0, 0, 0])
            _ = self.system.part.add(pos=[0, 0, 1])
            return p

        system = self.system

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
        system = self.system
        kT = 0.
        gamma = 1.42
        system.thermostat.set_dpd(kT=kT, seed=42)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.2,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.4)

        p0 = system.part.add(pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, .8, .3])
        p1 = system.part.add(pos=[3, 5, 5], type=0, v=v)

        system.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in system.part.all().f:
            np.testing.assert_array_equal(np.copy(f), [0., 0., 0.])

        # Only trans
        p1.pos = [5. - 1.3, 5, 5]

        system.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertLess(abs(p0.f[0]), 1e-14)
        np.testing.assert_allclose(
            np.copy(p0.f[1:2]), gamma * v[1:2], rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

        # Trans and parallel
        p1.pos = [5. - 1.1, 5, 5]

        system.integrator.run(0)

        np.testing.assert_allclose(
            np.copy(p0.f), gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

    def test_linear_weight_function(self):
        system = self.system
        kT = 0.
        gamma = 1.42
        system.thermostat.set_dpd(kT=kT, seed=42)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=1.2,
            trans_weight_function=1, trans_gamma=gamma, trans_r_cut=1.4)

        def calc_omega(dist, r_cut):
            return 1. - dist / r_cut

        p0 = system.part.add(pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, .8, .3])
        p1 = system.part.add(pos=[3, 5, 5], type=0, v=v)

        system.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in system.part.all().f:
            np.testing.assert_array_equal(np.copy(f), [0., 0., 0.])

        # Only trans
        p1.pos = [5. - 1.3, 5, 5]

        system.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertLess(abs(p0.f[0]), 1e-14)
        omega = calc_omega(1.3, 1.4)**2
        np.testing.assert_allclose(
            np.copy(p0.f[1:2]), omega * gamma * v[1:2], rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

        # Trans and parallel
        p1.pos = [5. - 1.1, 5, 5]

        system.integrator.run(0)

        omega = np.array([calc_omega(1.1, x)**2 for x in [1.2, 1.4, 1.4]])
        np.testing.assert_allclose(
            np.copy(p0.f), omega * gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

        # Trans and parallel 2nd point
        p1.pos = [5. - 0.5, 5, 5]

        system.integrator.run(0)

        omega = np.array([calc_omega(0.5, x)**2 for x in [1.2, 1.4, 1.4]])
        np.testing.assert_allclose(
            np.copy(p0.f), omega * gamma * v, rtol=0, atol=1e-11)
        np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

    def test_parabolic_weight_function(self):
        system = self.system
        kT = 0.
        gamma = 1.42
        kappa = 7.8
        r_cut = 1.2
        system.thermostat.set_dpd(kT=kT, seed=42)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, k=kappa, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=0.0, trans_r_cut=0.0)

        def calc_omega(dist, r_cut):
            return (1. - (dist / r_cut) ** kappa)

        p0 = system.part.add(pos=[5, 5, 5], type=0, v=[0, 0, 0])
        v = np.array([.5, 0., 0.])
        p1 = system.part.add(pos=[3, 5, 5], type=0, v=v)

        # Outside of both cutoffs, forces should be 0
        for f in system.part.all().f:
            np.testing.assert_array_equal(np.copy(f), [0., 0., 0.])

        # Place the particle at different positions to test the parabolic
        # weight function
        for dist in np.arange(0.1, 1.2, 50):

            p1.pos = [5. + dist, 5., 5.]
            system.integrator.run(0)
            omega = calc_omega(dist, r_cut)**2

            # The particle is moved along the x-direction. Hence, we are
            # testing the x element.
            np.testing.assert_allclose(
                np.copy(p0.f), omega * gamma * v, rtol=0, atol=1e-11)
            np.testing.assert_array_equal(np.copy(p0.f), -np.copy(p1.f))

    def test_ghosts_have_v(self):
        system = self.system

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
            system.part.add(pos=pos, v=v)

        gamma = 1.0
        system.thermostat.set_dpd(kT=0.0, seed=42)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=r_cut,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

        system.integrator.run(0)

        for p, ind in zip(system.part, ind_combinations):
            for i in ind:
                if ind[i] == 0:
                    sgn = 1
                else:
                    sgn = -1
                self.assertAlmostEqual(sgn * 4.0, p.f[i])

    def test_constraint(self):
        system = self.system

        dpd_vel = [1., 2., 3.]
        wall = espressomd.shapes.Wall(dist=1., normal=[1., 0., 0.])
        system.constraints.add(shape=wall, penetrable=True, particle_type=0,
                               particle_velocity=dpd_vel)

        system.thermostat.set_dpd(kT=0.0, seed=42)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=1., r_cut=1.0,
            trans_weight_function=0, trans_gamma=1., trans_r_cut=1.0)

        p1 = system.part.add(pos=[0.5, 0., 0.], type=0)
        p2 = system.part.add(pos=[1.5, 0., 0.], type=0)
        p3 = system.part.add(pos=[1.0, 0., 0.], type=0)
        system.integrator.run(0)
        np.testing.assert_array_almost_equal(np.copy(p1.f), dpd_vel)
        np.testing.assert_array_almost_equal(np.copy(p2.f), dpd_vel)
        np.testing.assert_array_almost_equal(np.copy(p3.f), 0.)

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

        system = self.system

        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=gamma / 2.0, trans_r_cut=r_cut)

        pos = system.box_l * np.random.random((n_part, 3))
        partcls = system.part.add(pos=pos)
        system.integrator.run(10)

        for kT in [0., 2.]:
            system.thermostat.set_dpd(kT=kT, seed=3)
            # run 1 integration step to get velocities
            partcls.v = np.zeros((n_part, 3))
            system.integrator.run(steps=1)

            pairs = system.part.pairs()

            stress = np.zeros([3, 3])

            for pair in pairs:
                dist = system.distance_vec(pair[0], pair[1])
                if np.linalg.norm(dist) < r_cut:
                    vel_diff = pair[1].v - pair[0].v
                    stress += calc_stress(dist, vel_diff)

            stress /= system.volume()

            dpd_stress = system.analysis.dpd_stress()

            dpd_obs = espressomd.observables.DPDStress()
            obs_stress = dpd_obs.calculate()

            np.testing.assert_array_almost_equal(np.copy(dpd_stress), stress)
            np.testing.assert_array_almost_equal(np.copy(obs_stress), stress)

    def test_momentum_conservation(self):
        r_cut = 1.0
        gamma = 5.
        r_cut = 2.9

        system = self.system
        system.thermostat.set_dpd(kT=1.3, seed=42)
        partcls = system.part.add(pos=((0, 0, 0), (0.1, 0.1, 0.1), (0.1, 0, 0)),
                                  mass=(1, 2, 3))

        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=r_cut,
            trans_weight_function=1, trans_gamma=gamma / 2.0, trans_r_cut=r_cut)
        momentum = np.matmul(partcls.v.T, partcls.mass)
        for _ in range(10):
            system.integrator.run(25)
            np.testing.assert_almost_equal(np.sum(partcls.f), 3 * [0.])
            np.testing.assert_allclose(
                np.matmul(partcls.v.T, partcls.mass), momentum, atol=1E-12)


if __name__ == "__main__":
    ut.main()
