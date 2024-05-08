#
# Copyright (C) 2010-2022 The ESPResSo project
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

import espressomd.lb
import espressomd.shapes
import espressomd.propagation
import unittest as ut
import unittest_decorators as utx
import numpy as np


class LBThermostatCommon:

    system = espressomd.System(box_l=[6., 6., 6.])
    system.time_step = 0.01
    system.cell_system.skin = 0.1
    system.min_global_cut = 2.

    def setUp(self):
        self.system.time_step = 0.01
        self.system.integrator.set_vv()

    def tearDown(self):
        self.system.lb = None
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def test_01__rng(self):
        """Test for RNG consistency."""

        system = self.system
        system.time_step = 0.01
        kT = 1.1
        gamma = 3.5

        def reset_fluid():
            lb_fluid = self.lb_class(agrid=1., tau=0.01, density=1., kT=kT,
                                     kinematic_viscosity=1., **self.lb_params)
            self.system.lb = lb_fluid
            return lb_fluid

        def reset_particle():
            self.system.part.clear()
            return system.part.add(pos=[0, 0, 0])

        self.assertIsNone(system.thermostat.kT)
        self.assertFalse(system.thermostat.lb.is_active)

        lb_fluid = reset_fluid()
        system.thermostat.set_lb(LB_fluid=lb_fluid, gamma=gamma, seed=41)
        system.thermostat.lb.call_method("override_philox_counter", counter=0)
        system.integrator.set_vv()

        self.assertIsNotNone(system.thermostat.kT)
        self.assertTrue(system.thermostat.lb.is_active)
        np.testing.assert_almost_equal(system.thermostat.kT, kT)
        np.testing.assert_almost_equal(
            np.copy(system.thermostat.lb.gamma), gamma)

        # run(0) does not increase the philox counter and should give the same
        # force
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force0 = np.copy(p.f)
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force1 = np.copy(p.f)
        np.testing.assert_almost_equal(force0, force1)
        np.testing.assert_almost_equal(force0, [0., 0., 0.])

        # run(1) should give a different force
        p = reset_particle()
        system.integrator.run(1)
        self.assertEqual(system.thermostat.lb.seed, 41)
        self.assertEqual(system.thermostat.lb.philox_counter, 1)
        force2 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force1, force2)))

        # Different seed should give a different force with same counter state
        # force2: lb.rng_counter() = 2, lb.rng_seed() = 41
        # force3: lb.rng_counter() = 2, lb.rng_seed() = 42
        reset_fluid()
        p = reset_particle()
        system.integrator.run(1)
        self.assertEqual(system.thermostat.lb.seed, 41)
        self.assertEqual(system.thermostat.lb.philox_counter, 2)
        force2 = np.copy(p.f)
        lb_fluid = reset_fluid()
        system.thermostat.set_lb(LB_fluid=lb_fluid, gamma=gamma, seed=42)
        system.thermostat.lb.call_method("override_philox_counter", counter=1)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.lb.seed, 42)
        self.assertEqual(system.thermostat.lb.philox_counter, 2)
        force3 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force2, force3)))

        # Same seed should not give the same force with different counter state
        # force3: lb.rng_counter() = 2, lb.rng_seed() = 42
        # force4: lb.rng_counter() = 3, lb.rng_seed() = 42
        lb_fluid = reset_fluid()
        p = reset_particle()
        system.thermostat.set_lb(LB_fluid=lb_fluid, gamma=gamma, seed=42)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.lb.seed, 42)
        self.assertEqual(system.thermostat.lb.philox_counter, 3)
        force4 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force3, force4)))

        # Seed offset should not give the same force with a lag
        # force4: lb.rng_counter() = 3, lb.rng_seed() = 42
        # force5: lb.rng_counter() = 4, lb.rng_seed() = 41
        lb_fluid = reset_fluid()
        reset_particle()
        system.thermostat.set_lb(
            LB_fluid=lb_fluid, kT=kT, gamma=gamma, seed=41)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.lb.seed, 41)
        self.assertEqual(system.thermostat.lb.philox_counter, 4)
        force5 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force4, force5)))

        self.system.thermostat.turn_off()
        self.assertFalse(system.thermostat.langevin.is_active)
        self.assertIsNone(system.thermostat.kT)

    @utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE",
                                "VIRTUAL_SITES_INERTIALESS_TRACERS"])
    def test_virtual_sites_relative(self):
        Propagation = espressomd.propagation.Propagation

        system = self.system
        lb_fluid = self.lb_class(agrid=1., tau=0.01, density=1.,
                                 kinematic_viscosity=1., **self.lb_params)
        system.lb = lb_fluid

        particle_inert = system.part.add(pos=[1, 1, 1], v=[0, 0, 0])
        particle_skip = system.part.add(pos=[1, 1, 1], v=[1, 0, 0],
                                        propagation=Propagation.NONE)
        tracer = system.part.add(pos=[0, 0, 0], v=[1, 0, 0],
                                 propagation=Propagation.TRANS_LB_TRACER)
        physical = system.part.add(pos=[0, 0, 0], v=[1, 0, 0])
        central_coupled = system.part.add(pos=[1, 1, 1], v=[1, 0, 0])
        virtual_coupled = system.part.add(pos=[2, 0, 0], v=[1, 0, 0])
        virtual_coupled.vs_auto_relate_to(central_coupled, couple_to_lb=True)
        virtual_decoupled = system.part.add(pos=[2, 0, 0], v=[5, 0, 0])
        virtual_decoupled.vs_auto_relate_to(
            central_coupled, couple_to_lb=False)
        central_recouple = system.part.add(pos=[1, 1, 1], v=[3, 0, 0],
                                           propagation=Propagation.NONE)
        virtual_coupler = system.part.add(pos=[2, 0, 0], v=[1, 0, 0])
        virtual_coupler.vs_auto_relate_to(central_recouple, couple_to_lb=True)

        system.thermostat.set_lb(LB_fluid=lb_fluid, gamma=1.)
        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(particle_inert.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(particle_skip.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(tracer.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])
        # vs relative: setup with additive forces
        np.testing.assert_almost_equal(np.copy(central_coupled.f), [-2, 0, 0])
        np.testing.assert_almost_equal(np.copy(virtual_coupled.f), [-1, 0, 0])
        np.testing.assert_almost_equal(np.copy(virtual_decoupled.f), [0, 0, 0])
        # vs relative: setup with overriden forces
        np.testing.assert_almost_equal(np.copy(central_recouple.f), [-3, 0, 0])
        np.testing.assert_almost_equal(np.copy(virtual_coupler.f), [-3, 0, 0])

    @utx.skipIfMissingFeatures(["MASS",
                                "ROTATIONAL_INERTIA",
                                "EXTERNAL_FORCES"])
    def test_propagation(self):
        Propagation = espressomd.propagation.Propagation
        np.random.seed(42)
        rtol_vel = 2e-4 if self.lb_params["single_precision"] else 1e-6
        v0 = np.array([-2., 2., 4.])
        o0 = np.array([2., 3., 4.])
        gamma_trans = 1.2
        gamma_rot = 0.6
        gamma_lbf = 4.
        positions = []

        def calc_trajectory(p, r0):
            t = self.system.time
            if p.propagation & (Propagation.SYSTEM_DEFAULT |
                                Propagation.TRANS_LB_MOMENTUM_EXCHANGE):
                friction_vel = np.exp(-gamma_lbf / p.mass * (t - tau / 2.))
                friction_pos = np.exp(-gamma_lbf / p.mass * t)
                x_drift = v0 * p.mass / gamma_lbf
                ref_vel = v0 * friction_vel
                ref_pos = r0 + x_drift * (1. - friction_pos)
            elif p.propagation & Propagation.TRANS_LANGEVIN:
                friction = np.exp(-gamma_trans / p.mass * t)
                x_drift = v0 * p.mass / gamma_trans
                ref_vel = v0 * friction
                ref_pos = r0 + x_drift * (1. - friction)
            elif p.propagation & Propagation.TRANS_NEWTON:
                ref_vel = v0
                ref_pos = r0 + v0 * t
            else:
                ref_vel = v0
                ref_pos = r0
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.ROT_LANGEVIN)):
                friction = np.exp(-gamma_rot / p.rinertia * t)
                o_term = p.ext_torque / gamma_rot
                ref_rot = o_term + (o0 - o_term) * friction
            elif p.propagation & Propagation.ROT_EULER:
                ref_rot = o0 + p.ext_torque / p.rinertia * t
            else:
                ref_rot = o0
            return np.copy(ref_pos), np.copy(ref_vel), np.copy(ref_rot)

        def make_particle(propagation):
            r0 = np.random.random(3) * self.system.box_l
            positions.append(r0)
            rotation = (propagation & (Propagation.SYSTEM_DEFAULT |
                        Propagation.ROT_EULER | Propagation.ROT_LANGEVIN)) != 0
            self.system.part.add(pos=r0, v=v0, propagation=propagation,
                                 omega_lab=o0, rotation=3 * [rotation],
                                 mass=4., rinertia=[1.5, 1.5, 1.5],
                                 ext_torque=[+0.25, -0.5, +1.])

        system = self.system
        system.time_step = tau = 0.00001
        system.thermostat.set_langevin(
            kT=0., gamma=gamma_trans, gamma_rotation=gamma_rot, seed=42)
        lb_fluid = self.lb_class(agrid=1., tau=tau, density=200.,
                                 kinematic_viscosity=1., **self.lb_params)
        system.lb = lb_fluid
        system.thermostat.set_lb(LB_fluid=lb_fluid, gamma=gamma_lbf)
        system.integrator.set_vv()
        make_particle(Propagation.NONE)
        make_particle(Propagation.SYSTEM_DEFAULT)
        make_particle(Propagation.TRANS_LANGEVIN | Propagation.ROT_LANGEVIN)
        make_particle(Propagation.TRANS_NEWTON | Propagation.ROT_EULER)
        make_particle(Propagation.TRANS_LB_MOMENTUM_EXCHANGE)
        make_particle(Propagation.TRANS_LB_MOMENTUM_EXCHANGE |
                      Propagation.ROT_EULER)
        make_particle(Propagation.TRANS_LB_MOMENTUM_EXCHANGE |
                      Propagation.ROT_LANGEVIN)
        system.time = 0.
        for power in range(self.lb_geom_progression):
            system.integrator.run(2**power)
            for p, r0 in zip(system.part.all(), positions):
                pos = np.copy(p.pos)
                vel = np.copy(p.v)
                rot = np.copy(p.omega_lab)
                ref_pos, ref_vel, ref_rot = calc_trajectory(p, r0)
                np.testing.assert_allclose(pos, ref_pos, rtol=1e-6)
                np.testing.assert_allclose(vel, ref_vel, rtol=rtol_vel)
                np.testing.assert_allclose(rot, ref_rot, rtol=1e-7)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaThermostatDoublePrecisionCPU(LBThermostatCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    lb_geom_progression = 10


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaThermostatSinglePrecisionCPU(LBThermostatCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    lb_geom_progression = 10


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBWalberlaThermostatDoublePrecisionGPU(LBThermostatCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": False}
    lb_geom_progression = 9


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBWalberlaThermostatSinglePrecisionGPU(LBThermostatCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": True}
    lb_geom_progression = 9


if __name__ == "__main__":
    ut.main()
