#
# Copyright (C) 2013-2022 The ESPResSo project
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
import espressomd
import espressomd.propagation
import itertools
import numpy as np


class BrownianThermostat(ut.TestCase):

    """Test Brownian Dynamics"""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.set_regular_decomposition(use_verlet_lists=True)
    system.cell_system.skin = 0
    system.periodicity = [False, False, False]

    def setUp(self):
        np.random.seed(42)

    def tearDown(self):
        self.system.time_step = 1e-12
        self.system.cell_system.skin = 0.0
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def check_rng(self, per_particle_gamma=False):
        """Test for RNG consistency."""

        kT = 1.1
        gamma = 3.5

        def reset_particle():
            self.system.part.clear()
            p = system.part.add(pos=[0, 0, 0])
            if espressomd.has_features("ROTATION"):
                p.rotation = [True, True, True]
            if per_particle_gamma:
                assert espressomd.has_features("THERMOSTAT_PER_PARTICLE")
                if espressomd.has_features("PARTICLE_ANISOTROPY"):
                    p.gamma = 3 * [gamma / 2]
                else:
                    p.gamma = gamma / 2
                if espressomd.has_features("ROTATION"):
                    p.gamma_rot = p.gamma * 1.5
            return p

        system = self.system
        system.time_step = 0.01

        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        system.integrator.set_brownian_dynamics()

        pos2force = np.sqrt(2 * kT / gamma * system.time_step)
        omega2torque = np.sqrt(2 * kT / gamma * system.time_step)

        # run(0) does not increase the philox counter and should give no force
        p = reset_particle()
        system.integrator.run(0)
        force0 = np.copy(p.pos) / pos2force
        np.testing.assert_almost_equal(force0, 0)

        # run(1) should give a force
        p = reset_particle()
        system.integrator.run(1)
        force1 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force1, [0, 0, 0])))
        if espressomd.has_features("ROTATION"):
            torque1 = np.copy(p.omega_body) / omega2torque
            self.assertTrue(np.all(np.not_equal(torque1, [0, 0, 0])))

        # Same seed should not give the same force with different counter state
        # force1: brownian.rng_counter() = 1, brownian.rng_seed() = 41
        # force2: brownian.rng_counter() = 2, brownian.rng_seed() = 41
        p = reset_particle()
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        system.integrator.run(1)
        force2 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force2, force1)))
        if espressomd.has_features("ROTATION"):
            torque2 = np.copy(p.omega_body) / omega2torque
            self.assertTrue(np.all(np.not_equal(torque2, torque1)))

        # Seed offset should not give the same force with a lag
        # force3: brownian.rng_counter() = 3, brownian.rng_seed() = 42
        # force4: brownian.rng_counter() = 4, brownian.rng_seed() = 41
        p = reset_particle()
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=42)
        system.integrator.run(1)
        force3 = np.copy(p.pos) / pos2force
        if espressomd.has_features("ROTATION"):
            torque3 = np.copy(p.omega_body) / omega2torque
        p = reset_particle()
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        system.integrator.run(1)
        force4 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force3, force4)))
        if espressomd.has_features("ROTATION"):
            torque4 = np.copy(p.omega_body) / omega2torque
            self.assertTrue(np.all(np.not_equal(torque3, torque4)))

    def test_01__rng(self):
        """Test for RNG consistency."""
        # No seed should throw exception
        with self.assertRaises(ValueError):
            self.system.thermostat.set_brownian(kT=1, gamma=2)
        self.check_rng()

    @utx.skipIfMissingFeatures("THERMOSTAT_PER_PARTICLE")
    def test_01__rng_per_particle(self):
        """Test for RNG consistency."""
        self.check_rng(True)

    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01
        Propagation = espressomd.propagation.Propagation

        virtual = system.part.add(pos=[0, 0, 0], v=[1, 0, 0],
                                  propagation=Propagation.NONE)
        physical = system.part.add(pos=[0, 0, 0], v=[1, 0, 0])

        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., seed=41)
        system.integrator.set_brownian_dynamics()

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [1, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])

    @utx.skipIfMissingFeatures(["ROTATION", "EXTERNAL_FORCES"])
    def test_fix_rotation(self):
        system = self.system
        system.time_step = 0.01
        part = system.part.add(
            pos=3 * [0.],
            fix=3 * [True],
            rotation=3 * [True])

        # torque only
        part.ext_torque = [0, 0, 1.3]
        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1.5, seed=41)
        system.integrator.set_brownian_dynamics()
        system.integrator.run(3)
        np.testing.assert_allclose(
            np.copy(part.omega_lab), [0, 0, 1.3 / 1.5], atol=1e-14)

        # noise only
        part.ext_torque = 3 * [0.]
        system.thermostat.set_brownian(
            kT=1, gamma=1, gamma_rotation=1.5, seed=41)
        system.integrator.run(3)
        self.assertGreater(np.linalg.norm(np.copy(part.omega_lab)), 0.)

    @utx.skipIfMissingFeatures(["MASS",
                                "ROTATIONAL_INERTIA",
                                "EXTERNAL_FORCES"])
    def test_propagation(self):
        """
        Check integration of Brownian's equations of motion and rotation
        for various combinations of propagation modes.
        """
        Propagation = espressomd.propagation.Propagation
        aniso = espressomd.has_features("PARTICLE_ANISOTROPY")
        gamma_trans = np.array([1.2, 1.2, 1.2]) if aniso else 1.2
        gamma_rot = np.array([1.5, 1.5, 1.5]) if aniso else 1.5
        x0 = np.array([0., 1., 2.])
        v0 = np.array([1., 2., 3.])
        o0 = np.array([2., 3., 4.])
        ext_force = np.array([-1., +2., -4.])
        ext_torque = np.array([+1., -3., +5.])

        def has_gamma(gamma):
            return gamma[0] > 0. if np.shape(gamma) else gamma > 0.

        def calc_trajectory(p):
            t = self.system.time
            gamma_t = gamma_trans
            gamma_r = gamma_rot
            if espressomd.has_features("THERMOSTAT_PER_PARTICLE") and \
                    p.propagation & Propagation.SYSTEM_DEFAULT:
                if has_gamma(p.gamma):
                    gamma_t = p.gamma
                if has_gamma(p.gamma_rot):
                    gamma_r = p.gamma_rot
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.TRANS_BROWNIAN)):
                ref_vel = p.ext_force / gamma_t
                ref_pos = x0 + p.ext_force / gamma_t * t
            elif p.propagation & Propagation.TRANS_NEWTON:
                ref_vel = v0 + p.ext_force / p.mass * t
                ref_pos = x0 + v0 * t + 0.5 * p.ext_force / p.mass * t**2
            else:
                ref_vel = v0
                ref_pos = x0
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.ROT_BROWNIAN)):
                ref_rot = p.ext_torque / gamma_r
            elif p.propagation & Propagation.ROT_EULER:
                ref_rot = o0 + p.ext_torque / p.rinertia * t
            else:
                ref_rot = o0
            return np.copy(ref_pos), np.copy(ref_vel), np.copy(ref_rot)

        system = self.system
        system.time_step = 0.0001
        system.thermostat.set_brownian(
            kT=0., gamma=gamma_trans, gamma_rotation=gamma_rot, seed=42)
        system.integrator.set_brownian_dynamics()
        modes_trans = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.TRANS_NEWTON,
            Propagation.TRANS_BROWNIAN,
        ]
        modes_rot = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.ROT_EULER,
            Propagation.ROT_BROWNIAN,
        ]
        for mode_trans, mode_rot in itertools.product(modes_trans, modes_rot):
            try:
                rotation = (mode_rot != Propagation.NONE) or \
                           (mode_trans == Propagation.SYSTEM_DEFAULT)
                system.part.add(pos=x0, v=v0, omega_lab=o0,
                                ext_force=ext_force, ext_torque=ext_torque,
                                rotation=3 * [rotation],
                                mass=0.5, rinertia=[0.5, 0.5, 0.5],
                                propagation=mode_trans | mode_rot)
            except BaseException:
                continue
        assert len(system.part) > 3
        if espressomd.has_features("THERMOSTAT_PER_PARTICLE"):
            ps = system.part.select(propagation=Propagation.SYSTEM_DEFAULT)
            list(ps)[1].gamma = 0.8 * gamma_trans
            list(ps)[2].gamma_rot = 0.7 * gamma_rot
        system.time = 0.
        for i in range(10):
            system.integrator.run(2**i)
            for p in system.part.all():
                pos = np.copy(p.pos)
                vel = np.copy(p.v)
                rot = np.copy(p.omega_lab)
                ref_pos, ref_vel, ref_rot = calc_trajectory(p)
                np.testing.assert_allclose(pos, ref_pos, rtol=1e-10)
                np.testing.assert_allclose(vel, ref_vel, rtol=1e-10)
                np.testing.assert_allclose(rot, ref_rot, rtol=1e-10)


if __name__ == "__main__":
    ut.main()
