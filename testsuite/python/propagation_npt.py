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
import tests_common
import numpy as np
import espressomd
import espressomd.integrate
import espressomd.propagation


@utx.skipIfMissingFeatures(["NPT"])
class IntegratorNPT(ut.TestCase):

    """This tests the NpT integrator interface."""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.box_l = [5] * 3
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.25

    def tearDown(self):
        self.system.part.clear()
        if espressomd.has_features(["ELECTROSTATICS"]):
            self.system.electrostatics.clear()
        if espressomd.has_features(["DIPOLES"]):
            self.system.magnetostatics.clear()
        # reset RNG counter to make tests independent of execution order
        self.system.thermostat.npt_iso.call_method(
            "override_philox_counter", counter=0)
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        if espressomd.has_features(["LENNARD_JONES"]):
            self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()

    def test_integrator_exceptions(self):
        # invalid parameters should throw exceptions
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=-1, piston=1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=1, piston=0)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=1, piston=-1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1, piston=1, direction=[False, False, False])
        with self.assertRaises(Exception):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1, piston=1, direction=[True, False])

    @utx.skipIfMissingFeatures(["MASS", "EXTERNAL_FORCES"])
    def test_00_propagation(self):
        """
        Check integration of the equations of motion.
        """
        Propagation = espressomd.propagation.Propagation
        gamma0 = 1.2
        v0 = np.array([1., 2., 3.])
        ext_force = np.array([-2., +2., -3.])

        def calc_trajectory(p, x0):
            t = self.system.time
            gamma_t = gamma0
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.TRANS_LANGEVIN_NPT)):
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
            return np.copy(ref_pos), np.copy(ref_vel)

        system = self.system
        system.time_step = 0.00001
        system.thermostat.set_npt(kT=0., gamma0=gamma0, gammav=1e-6, seed=42)
        system.integrator.set_isotropic_npt(ext_pressure=0.01, piston=1e6)
        positions = []
        modes_trans = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.TRANS_NEWTON,
            Propagation.TRANS_LANGEVIN_NPT,
        ]
        for mode_trans in modes_trans:
            x0 = np.random.random(3) * np.copy(system.box_l)
            system.part.add(pos=x0, v=v0, ext_force=ext_force,
                            mass=4., propagation=mode_trans)
            positions.append(x0)
        system.time = 0.
        for i in range(10):
            system.integrator.run(2**i)
            for p, x0 in zip(system.part.all(), positions):
                pos = np.copy(p.pos)
                vel = np.copy(p.v)
                ref_pos, ref_vel = calc_trajectory(p, x0)
                np.testing.assert_allclose(pos, ref_pos, rtol=1e-7)
                np.testing.assert_allclose(vel, ref_vel, rtol=1e-7)

    @utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE"])
    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01

        virtual = system.part.add(pos=[0, 0, 0], v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], v=[2, 0, 0])
        virtual.vs_relative = (physical.id, 0.1, (1., 0., 0., 0.))

        system.thermostat.set_npt(kT=0., gamma0=2., gammav=1e-6, seed=42)
        system.integrator.set_isotropic_npt(ext_pressure=0.01, piston=1e6)

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(physical.f), [0., 0., 0.])
        np.testing.assert_almost_equal(np.copy(physical.v), [1.9602, 0., 0.])
        np.testing.assert_almost_equal(np.copy(physical.pos), [0.0198, 0., 0.])
        np.testing.assert_almost_equal(np.copy(virtual.f), [0., 0., 0.])
        np.testing.assert_almost_equal(np.copy(virtual.v), [1.9602, 0., 0.])
        np.testing.assert_almost_equal(np.copy(virtual.pos), [0.0198, 0., 0.1])

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_09_integrator_recovery(self):
        # the system is still in a valid state after a failure
        system = self.system
        np.random.seed(42)
        npt_params = {'ext_pressure': 0.01, 'piston': 0.001}
        system.box_l = [6] * 3
        system.part.add(pos=np.random.uniform(0, system.box_l[0], (11, 3)))
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1, cutoff=2**(1 / 6), shift=0.25)
        system.thermostat.set_npt(kT=1.0, gamma0=2, gammav=0.04, seed=42)
        system.integrator.set_isotropic_npt(**npt_params)

        # get the equilibrium box length for the chosen NpT parameters
        system.integrator.run(500)
        # catch unstable simulation early (when the DP3M test case ran first)
        assert system.box_l[0] < 20., "NpT simulation is unstable"
        system.integrator.run(1500)
        box_l_ref = system.box_l[0]

        # resetting the NpT integrator with incorrect values doesn't leave the
        # system in an undefined state (the old parameters aren't overwritten)
        with self.assertRaises(RuntimeError):
            system.integrator.set_isotropic_npt(ext_pressure=-1, piston=100)
        with self.assertRaises(RuntimeError):
            system.integrator.set_isotropic_npt(ext_pressure=100, piston=-1)
        # the core state is unchanged
        system.integrator.run(500)
        self.assertAlmostEqual(system.box_l[0], box_l_ref, delta=0.15)

        # setting another integrator with incorrect values doesn't leave the
        # system in an undefined state (the old integrator is still active)
        with self.assertRaises(RuntimeError):
            system.integrator.set_steepest_descent(
                f_max=-10, gamma=0, max_displacement=0.1)
        # the interface state is unchanged
        integrator_state = system.integrator.get_params()
        self.assertIsInstance(integrator_state['integrator'],
                              espressomd.integrate.VelocityVerletIsotropicNPT)
        params = integrator_state['integrator'].get_params()
        self.assertEqual(params['ext_pressure'], npt_params['ext_pressure'])
        self.assertEqual(params['piston'], npt_params['piston'])
        # the core state is unchanged
        system.integrator.run(500)
        self.assertAlmostEqual(system.box_l[0], box_l_ref, delta=0.15)

        # setting the NpT integrator with incorrect values doesn't leave the
        # system in an undefined state (the old integrator is still active)
        system.thermostat.turn_off()
        system.integrator.set_vv()
        system.part.clear()
        system.box_l = [5] * 3
        positions_start = np.array([[0, 0, 0], [1., 0, 0]])
        system.part.add(pos=positions_start)
        with self.assertRaises(RuntimeError):
            system.integrator.set_isotropic_npt(ext_pressure=-1, piston=100)
        # the interface state is unchanged
        self.assertIsInstance(system.integrator.get_params()['integrator'],
                              espressomd.integrate.VelocityVerlet)
        # the core state is unchanged
        system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(system.part.all().pos),
            positions_start + np.array([[-1.2e-3, 0, 0], [1.2e-3, 0, 0]]))

    def run_with_p3m(self, container, p3m, method):
        system = self.system
        npt_kwargs = {"ext_pressure": 0.001, "piston": 0.001}
        npt_kwargs_rectangular = {
            "cubic_box": False, "direction": (False, True, True), **npt_kwargs}
        np.random.seed(42)
        # set up particles
        system.box_l = [6] * 3
        partcl = system.part.add(
            pos=np.random.uniform(0, system.box_l[0], (11, 3)))
        if espressomd.has_features("P3M"):
            partcl.q = np.sign(np.arange(-5, 6))
        if espressomd.has_features("DP3M"):
            partcl.dip = tests_common.random_dipoles(11)
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1, cutoff=2**(1 / 6), shift=0.25)
        system.integrator.set_steepest_descent(
            f_max=10, gamma=0.1, max_displacement=0.01)
        system.integrator.run(100)
        system.integrator.set_vv()
        # combine NpT with a P3M algorithm
        container.solver = p3m
        system.integrator.run(10)
        system.integrator.set_isotropic_npt(**npt_kwargs)
        system.thermostat.set_npt(kT=1.0, gamma0=2, gammav=0.04, seed=42)
        system.integrator.run(10)
        # check runtime warnings
        system.thermostat.turn_off()
        system.integrator.set_vv()
        err_msg = f"If {method} is being used you must use the cubic box NpT"
        with self.assertRaisesRegex(RuntimeError, err_msg):
            system.integrator.set_isotropic_npt(**npt_kwargs_rectangular)
        self.assertIsInstance(
            system.integrator.integrator,
            espressomd.integrate.VelocityVerlet)
        container.solver = None
        system.thermostat.set_npt(kT=1.0, gamma0=2, gammav=0.04, seed=42)
        system.integrator.set_isotropic_npt(**npt_kwargs_rectangular)
        container.solver = p3m
        with self.assertRaisesRegex(Exception, err_msg):
            system.integrator.run(0, recalc_forces=True)

    @utx.skipIfMissingFeatures(["DP3M", "LENNARD_JONES"])
    def test_npt_dp3m_cpu(self):
        import espressomd.magnetostatics
        dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=1.0, accuracy=1e-2, mesh=3 * [36], cao=7, r_cut=1.0,
            alpha=2.995, tune=False)
        self.run_with_p3m(self.system.magnetostatics, dp3m, "magnetostatics")

    @utx.skipIfMissingFeatures(["P3M", "LENNARD_JONES"])
    def test_npt_p3m_cpu(self):
        import espressomd.electrostatics
        p3m = espressomd.electrostatics.P3M(
            prefactor=1.0, accuracy=1e-2, mesh=3 * [8], cao=3, r_cut=0.36,
            alpha=5.35, tune=False)
        self.run_with_p3m(self.system.electrostatics, p3m, "electrostatics")

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M", "LENNARD_JONES"])
    def test_npt_p3m_gpu(self):
        import espressomd.electrostatics
        p3m = espressomd.electrostatics.P3MGPU(
            prefactor=1.0, accuracy=1e-2, mesh=3 * [8], cao=3, r_cut=0.36,
            alpha=5.35, tune=False)
        self.run_with_p3m(self.system.electrostatics, p3m, "electrostatics")


if __name__ == "__main__":
    ut.main()
