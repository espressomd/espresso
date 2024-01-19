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
import espressomd.lb
import espressomd.propagation
import itertools
import numpy as np


class LangevinThermostat(ut.TestCase):

    """Test Langevin Dynamics"""
    system = espressomd.System(box_l=[12., 12., 12.])
    system.cell_system.set_regular_decomposition(use_verlet_lists=True)
    system.cell_system.skin = 0.
    system.min_global_cut = 2.
    system.periodicity = [False, False, False]

    def setUp(self):
        np.random.seed(42)
        self.system.time_step = 1e-12
        self.system.cell_system.skin = 0.0
        self.system.integrator.set_vv()

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        if espressomd.has_features("WALBERLA"):
            self.system.lb = None

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

        self.assertIsNone(system.thermostat.kT)
        self.assertFalse(system.thermostat.langevin.is_active)

        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=41)
        system.thermostat.langevin.call_method(
            "override_philox_counter", counter=0)
        system.integrator.set_vv()

        self.assertIsNotNone(system.thermostat.kT)
        np.testing.assert_almost_equal(system.thermostat.kT, kT)
        np.testing.assert_almost_equal(
            np.copy(system.thermostat.langevin.gamma), gamma)
        if espressomd.has_features("ROTATION"):
            np.testing.assert_almost_equal(
                np.copy(system.thermostat.langevin.gamma_rotation), gamma)

        # run(0) does not increase the philox counter and should give the same
        # force
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force0 = np.copy(p.f)
        if espressomd.has_features("ROTATION"):
            torque0 = np.copy(p.torque_lab)
        system.integrator.run(0, recalc_forces=True)
        force1 = np.copy(p.f)
        np.testing.assert_almost_equal(force0, force1)
        if espressomd.has_features("ROTATION"):
            torque1 = np.copy(p.torque_lab)
            np.testing.assert_almost_equal(torque0, torque1)

        # run(1) should give a different force
        p = reset_particle()
        system.integrator.run(1)
        self.assertEqual(system.thermostat.langevin.philox_counter, 1)
        force2 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force1, force2)))
        if espressomd.has_features("ROTATION"):
            torque2 = np.copy(p.torque_lab)
            self.assertTrue(np.all(np.not_equal(torque1, torque2)))

        # Different seed should give a different force with same counter state
        # force2: langevin.rng_counter() = 1, langevin.rng_seed() = 41
        # force3: langevin.rng_counter() = 1, langevin.rng_seed() = 42
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        self.assertEqual(system.thermostat.langevin.seed, 41)
        self.assertEqual(system.thermostat.langevin.philox_counter, 1)
        force2 = np.copy(p.f)
        if espressomd.has_features("ROTATION"):
            torque2 = np.copy(p.torque_lab)
        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
        system.integrator.run(0, recalc_forces=True)
        self.assertEqual(system.thermostat.langevin.seed, 42)
        self.assertEqual(system.thermostat.langevin.philox_counter, 1)
        force3 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force2, force3)))
        if espressomd.has_features("ROTATION"):
            torque3 = np.copy(p.torque_lab)
            self.assertTrue(np.all(np.not_equal(torque2, torque3)))

        # Same seed should not give the same force with different counter state
        p = reset_particle()
        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.langevin.seed, 42)
        self.assertEqual(system.thermostat.langevin.philox_counter, 2)
        force4 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force3, force4)))
        if espressomd.has_features("ROTATION"):
            torque4 = np.copy(p.torque_lab)
            self.assertTrue(np.all(np.not_equal(torque3, torque4)))

        # Seed offset should not give the same force with a lag
        # force4: langevin.rng_counter() = 2, langevin.rng_seed() = 42
        # force5: langevin.rng_counter() = 3, langevin.rng_seed() = 41
        reset_particle()
        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=41)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.langevin.seed, 41)
        self.assertEqual(system.thermostat.langevin.philox_counter, 3)
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force5 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force4, force5)))
        if espressomd.has_features("ROTATION"):
            torque5 = np.copy(p.torque_lab)
            self.assertTrue(np.all(np.not_equal(torque4, torque5)))

        with self.assertRaises(ValueError):
            system.thermostat.set_langevin(kT=-1., gamma=2.)
        with self.assertRaises(ValueError):
            system.thermostat.set_langevin(kT=1., gamma=-2.)

        self.system.thermostat.turn_off()
        self.assertFalse(system.thermostat.langevin.is_active)
        self.assertIsNone(system.thermostat.kT)

    def test_01__rng(self):
        """Test for RNG consistency."""
        # No seed should throw exception
        with self.assertRaises(ValueError):
            self.system.thermostat.set_langevin(kT=1, gamma=2)
        self.check_rng()

    @utx.skipIfMissingFeatures("THERMOSTAT_PER_PARTICLE")
    def test_01__rng_per_particle(self):
        """Test for RNG consistency."""
        self.check_rng(True)

    def test_02__friction_trans(self):
        """Test the translational friction-only part of the thermostat."""

        system = self.system
        v0 = np.array([5., 5., 5.])
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            gamma = np.array([0.5, 2., 1.5])
            atol = 3e-5
        else:
            gamma = 2.
            atol = 3e-4

        system.time = 0.
        system.time_step = 0.0001
        p = system.part.add(pos=(0, 0, 0), v=v0)
        if espressomd.has_features("MASS"):
            p.mass = 3
        system.thermostat.set_langevin(kT=0, gamma=gamma, seed=41)

        system.time = 0
        for _ in range(100):
            system.integrator.run(10)
            v_ref = v0 * np.exp(-gamma / p.mass * system.time)
            np.testing.assert_allclose(np.copy(p.v), v_ref, atol=atol)

    @utx.skipIfMissingFeatures("ROTATION")
    def test_03__friction_rot(self):
        """Test the rotational friction-only part of the thermostat."""

        system = self.system
        aniso = espressomd.has_features("PARTICLE_ANISOTROPY")
        o0 = np.array([5., 5., 5.])
        gamma_t = np.array([0.5, 2., 1.5]) if aniso else 2.
        gamma_r = np.array([1.5, 0.7, 1.2]) if aniso else 3.
        rinertia = np.array([1., 1., 1.])

        system.time = 0.
        system.time_step = 0.0001
        p = system.part.add(pos=(0, 0, 0), omega_body=o0, rotation=3 * [True])
        if espressomd.has_features("ROTATIONAL_INERTIA"):
            rinertia = np.array([2., 2., 2.])
            p.rinertia = rinertia
        system.thermostat.set_langevin(
            kT=0, gamma=gamma_t, gamma_rotation=gamma_r, seed=41)

        system.time = 0
        for _ in range(100):
            system.integrator.run(10)
            ref_omega_body = o0 * np.exp(-gamma_r / rinertia * system.time)
            np.testing.assert_allclose(
                np.copy(p.omega_body), ref_omega_body, atol=5E-4)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_07__virtual(self):
        system = self.system
        system.time_step = dt = 0.03

        v0 = np.array([2, 0, 0])
        virtual = system.part.add(pos=[0, 0, 0], v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], v=v0)
        virtual.vs_relative = (physical.id, 0.01, (1., 0., 0., 0.))

        system.thermostat.set_langevin(
            kT=0, gamma=1, gamma_rotation=1., seed=41)

        system.integrator.run(0)

        np.testing.assert_almost_equal(np.copy(virtual.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), -v0)

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), dt * v0 / 2. - v0)

    @utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE", "WALBERLA"])
    def test_virtual_sites_relative(self):
        Propagation = espressomd.propagation.Propagation
        system = self.system
        system.time_step = 0.01
        o0 = np.array([0, 0, 0])

        system.lb = espressomd.lb.LBFluidWalberla(
            tau=0.01, agrid=2., density=1., kinematic_viscosity=1., kT=0.)
        system.thermostat.set_lb(LB_fluid=system.lb, seed=42, gamma=1.)
        system.thermostat.set_langevin(
            kT=0., gamma=1., gamma_rotation=1., seed=42)

        real = system.part.add(pos=[1, 1, 1], v=[3, 4, 5], omega_lab=o0)
        virt_lb = system.part.add(pos=[2, 0, 0], v=[1, 2, 3], omega_lab=o0)
        virt_lg = system.part.add(pos=[2, 0, 0], v=[4, 5, 6], omega_lab=o0)
        virt_lx = system.part.add(pos=[2, 0, 0], v=[7, 8, 9], omega_lab=o0)
        system.part.all().propagation = (
            Propagation.TRANS_LANGEVIN | Propagation.ROT_LANGEVIN)

        refs = {
            "real.f": [-3, -4, -5], "real.torque_lab": [0, 0, 0],
            "virt_lb.f": [-1, -2, -3], "virt_lb.torque_lab": [0, 0, 0],
            "virt_lg.f": [-4, -5, -6], "virt_lg.torque_lab": [0, 0, 0],
            "virt_lx.f": [-7, -8, -9], "virt_lx.torque_lab": [0, 0, 0],
        }

        def check():
            assertion = np.testing.assert_almost_equal
            system.integrator.run(0, recalc_forces=True)
            assertion(np.copy(real.f), refs["real.f"])
            assertion(np.copy(virt_lg.f), refs["virt_lg.f"])
            assertion(np.copy(virt_lb.f), refs["virt_lb.f"])
            assertion(np.copy(virt_lx.f), refs["virt_lx.f"])
            assertion(np.copy(real.torque_lab), refs["real.torque_lab"])
            assertion(np.copy(virt_lg.torque_lab), refs["virt_lg.torque_lab"])
            assertion(np.copy(virt_lb.torque_lab), refs["virt_lb.torque_lab"])
            assertion(np.copy(virt_lx.torque_lab), refs["virt_lx.torque_lab"])

        check()

        virt_lg.vs_auto_relate_to(real, couple_to_langevin=True)
        self.assertEqual(virt_lg.propagation,
                         Propagation.TRANS_VS_RELATIVE |
                         Propagation.TRANS_LANGEVIN |
                         Propagation.ROT_VS_RELATIVE |
                         Propagation.ROT_LANGEVIN)
        refs["real.f"] = [-6, -8, -10]
        refs["real.torque_lab"] = [1, 8, -7]
        refs["virt_lg.f"] = [-3, -4, -5]
        check()

        virt_lb.vs_auto_relate_to(real, couple_to_lb=True)
        self.assertEqual(virt_lb.propagation,
                         Propagation.TRANS_LB_MOMENTUM_EXCHANGE |
                         Propagation.TRANS_VS_RELATIVE |
                         Propagation.ROT_VS_RELATIVE)
        refs["virt_lb.f"] = [0, 0, 0]
        check()

        virt_lx.vs_auto_relate_to(
            real, couple_to_lb=True, couple_to_langevin=True)
        self.assertEqual(virt_lx.propagation,
                         Propagation.TRANS_LB_MOMENTUM_EXCHANGE |
                         Propagation.TRANS_VS_RELATIVE |
                         Propagation.ROT_VS_RELATIVE |
                         Propagation.ROT_LANGEVIN)
        refs["virt_lx.f"] = [0, 0, 0]
        check()

        virt_lb.omega_lab = [1, 0, 0]
        refs["virt_lb.torque_lab"] = [0, 0, 0]
        refs["real.torque_lab"] = [1, 8, -7]
        check()

        virt_lg.omega_lab = [1, 0, 0]
        refs["virt_lg.torque_lab"] = [-1, 0, 0]
        refs["real.torque_lab"] = [0, 8, -7]
        check()

        virt_lx.omega_lab = [1, 0, 0]
        refs["virt_lx.torque_lab"] = [-1, 0, 0]
        refs["real.torque_lab"] = [-1, 8, -7]
        check()

        real.omega_lab = [1, 0, 0]
        refs["real.f"] = [-6, -9, -9]
        refs["real.torque_lab"] = [-4, 7, -8]
        refs["virt_lg.f"] = [-3, -5, -4]
        check()

    @utx.skipIfMissingFeatures(["MASS",
                                "ROTATIONAL_INERTIA",
                                "EXTERNAL_FORCES"])
    def test_propagation(self):
        """
        Check integration of Langevin's equations of motion and rotation
        for various combinations of propagation modes.
        """
        Propagation = espressomd.propagation.Propagation
        aniso = espressomd.has_features("PARTICLE_ANISOTROPY")
        gamma_trans = np.array([1.2, 1.2, 1.2]) if aniso else 1.2
        gamma_rot = np.array([0.6, 0.6, 0.6]) if aniso else 0.6
        x0 = np.array([0., 1., 2.])
        v0 = np.array([1., 2., 3.])
        o0 = np.array([1., 2., 1.])
        ext_force = np.array([-2., +2., -3.])
        ext_torque = np.array([+0.25, -0.5, +1.])

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
                                 Propagation.TRANS_LANGEVIN)):
                friction = np.exp(-gamma_t / p.mass * t)
                v_term = p.ext_force / gamma_t
                x_drift = (v0 - v_term) * p.mass / gamma_t
                ref_vel = v_term + (v0 - v_term) * friction
                ref_pos = x0 + v_term * t + x_drift * (1. - friction)
            elif p.propagation & Propagation.TRANS_NEWTON:
                ref_vel = v0 + p.ext_force / p.mass * t
                ref_pos = x0 + v0 * t + 0.5 * p.ext_force / p.mass * t**2
            else:
                ref_vel = v0
                ref_pos = x0
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.ROT_LANGEVIN)):
                friction = np.exp(-gamma_r / p.rinertia * t)
                o_term = p.ext_torque / gamma_r
                ref_rot = o_term + (o0 - o_term) * friction
            elif p.propagation & Propagation.ROT_EULER:
                ref_rot = o0 + p.ext_torque / p.rinertia * t
            else:
                ref_rot = o0
            return np.copy(ref_pos), np.copy(ref_vel), np.copy(ref_rot)

        system = self.system
        system.time_step = 0.00001
        system.thermostat.set_langevin(
            kT=0., gamma=gamma_trans, gamma_rotation=gamma_rot, seed=42)
        system.integrator.set_vv()
        modes_trans = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.TRANS_NEWTON,
            Propagation.TRANS_LANGEVIN,
        ]
        modes_rot = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.ROT_EULER,
            Propagation.ROT_LANGEVIN,
        ]
        for mode_trans, mode_rot in itertools.product(modes_trans, modes_rot):
            try:
                rotation = (mode_rot != Propagation.NONE) or \
                           (mode_trans == Propagation.SYSTEM_DEFAULT)
                system.part.add(pos=x0, v=v0, omega_lab=o0,
                                ext_force=ext_force, ext_torque=ext_torque,
                                rotation=3 * [rotation],
                                mass=4., rinertia=[1.5, 1.5, 1.5],
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
                np.testing.assert_allclose(pos, ref_pos, rtol=1e-7)
                np.testing.assert_allclose(vel, ref_vel, rtol=1e-7)
                np.testing.assert_allclose(rot, ref_rot, rtol=1e-7)


if __name__ == "__main__":
    ut.main()
