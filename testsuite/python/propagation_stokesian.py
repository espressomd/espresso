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


@utx.skipIfMissingFeatures("STOKESIAN_DYNAMICS")
class StokesianThermostat(ut.TestCase):

    """Test Stokesian Dynamics"""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.skin = 0
    system.periodicity = [False, False, False]

    def setUp(self):
        np.random.seed(42)

    def tearDown(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.0
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        self.system.periodicity = [False, False, False]

    def test_01__rng(self):
        """Test for RNG consistency."""
        def reset_particle():
            self.system.part.clear()
            p = self.system.part.add(pos=[0, 0, 0])
            return p

        system = self.system
        system.time_step = 0.01

        kT = 1e-4
        radius = 1.5
        viscosity = 2.4

        # invalid parameters should throw exceptions
        with self.assertRaisesRegex(ValueError, "Parameter 'kT' cannot be negative"):
            system.thermostat.set_stokesian(kT=-1)
        with self.assertRaises(ValueError):
            system.thermostat.set_stokesian(kT=1, seed=-1)

        system.thermostat.set_stokesian(kT=kT, seed=41)
        system.integrator.set_stokesian_dynamics(
            viscosity=viscosity, radii={0: radius})

        pos2force = np.sqrt(24 * kT * system.time_step)

        # run(0) does not increase the philox counter and should give no force
        p = reset_particle()
        system.integrator.run(0)
        self.assertEqual(system.thermostat.stokesian.seed, 41)
        self.assertEqual(system.thermostat.stokesian.philox_counter, 0)
        force0 = np.copy(p.pos) / pos2force
        np.testing.assert_almost_equal(force0, 0)

        # run(1) should give a force
        p = reset_particle()
        system.integrator.run(1)
        self.assertEqual(system.thermostat.stokesian.seed, 41)
        self.assertEqual(system.thermostat.stokesian.philox_counter, 1)
        force1 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force1, [0, 0, 0])))

        # Same seed should not give the same force with different counter state
        # force1: brownian.rng_counter() = 1, brownian.rng_seed() = 41
        # force2: brownian.rng_counter() = 2, brownian.rng_seed() = 41
        p = reset_particle()
        system.thermostat.set_stokesian(kT=kT, seed=41)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.stokesian.seed, 41)
        self.assertEqual(system.thermostat.stokesian.philox_counter, 2)
        force2 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force2, force1)))

        # Seed offset should not give the same force with a lag
        # force3: brownian.rng_counter() = 3, brownian.rng_seed() = 42
        # force4: brownian.rng_counter() = 4, brownian.rng_seed() = 41
        p = reset_particle()
        system.thermostat.set_stokesian(kT=kT, seed=42)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.stokesian.seed, 42)
        self.assertEqual(system.thermostat.stokesian.philox_counter, 3)
        force3 = np.copy(p.pos) / pos2force
        p = reset_particle()
        system.thermostat.set_stokesian(kT=kT, seed=41)
        system.integrator.run(1)
        self.assertEqual(system.thermostat.stokesian.seed, 41)
        self.assertEqual(system.thermostat.stokesian.philox_counter, 4)
        force4 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force3, force4)))

    def test_integrator_exceptions(self):
        # invalid parameters should throw exceptions
        with self.assertRaisesRegex(ValueError, "Particle radius for type 0 has an invalid value"):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=1.0, radii={0: -1})
        with self.assertRaisesRegex(ValueError, "Viscosity has an invalid value"):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=-1.0, radii={0: 1.0})
        with self.assertRaisesRegex(ValueError, "Unknown approximation 'STS'"):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=1.0, radii={0: 1.0}, approximation_method="STS")

        # invalid PBC should throw exceptions
        self.system.integrator.set_vv()
        self.system.periodicity = [False, False, True]
        with self.assertRaisesRegex(RuntimeError, r"Stokesian Dynamics requires periodicity \(False, False, False\)"):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=1.0, radii={0: 1.0})

        self.system.periodicity = [False, False, False]
        self.system.integrator.set_stokesian_dynamics(
            viscosity=1.0, radii={0: 1.0})

        with self.assertRaisesRegex(Exception, r"Stokesian Dynamics requires periodicity \(False, False, False\)"):
            self.system.periodicity = [False, True, False]

    @utx.skipIfMissingFeatures(["MASS",
                                "ROTATIONAL_INERTIA",
                                "EXTERNAL_FORCES"])
    def test_propagation(self):
        """
        Check integration of the equations of motion and rotation
        for various combinations of propagation modes.
        The exact values of Stokesian trajectories are not tested; instead
        the Stokesian thermostat parameters are chosen such that positions,
        velocities and angular velocities are in the range [-100, +100]
        in simulation units.
        """
        Propagation = espressomd.propagation.Propagation
        np.random.seed(42)
        x0 = [0., 1., 2.]
        v0 = np.array([1., 2., 3.])
        o0 = np.array([2., 3., 4.])
        ext_force = np.array([-1., +2., -4.])
        ext_torque = np.array([+1., -3., +5.])

        def calc_trajectory(p):
            t = self.system.time
            atol = 1e-6
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.TRANS_STOKESIAN)):
                atol = 1e3
                ref_pos = np.array([0., 0., 0.])
                ref_vel = np.array([0., 0., 0.])
            elif p.propagation & Propagation.TRANS_NEWTON:
                ref_vel = v0 + p.ext_force / p.mass * t
                ref_pos = x0 + v0 * t + 0.5 * p.ext_force / p.mass * t**2
            else:
                ref_vel = v0
                ref_pos = x0
            if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                 Propagation.ROT_STOKESIAN)):
                atol = 1e3
                ref_rot = np.array([0., 0., 0.])
            elif p.propagation & Propagation.ROT_EULER:
                ref_rot = o0 + p.ext_torque / p.rinertia * t
            else:
                ref_rot = o0
            return np.copy(ref_pos), np.copy(ref_vel), np.copy(ref_rot), atol

        system = self.system
        system.time_step = 0.0001
        system.thermostat.set_stokesian(kT=1., seed=42)
        system.integrator.set_stokesian_dynamics(viscosity=2.4, radii={0: 1.5})
        modes_trans = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.TRANS_NEWTON,
            Propagation.TRANS_STOKESIAN,
        ]
        modes_rot = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.ROT_EULER,
            Propagation.ROT_STOKESIAN,
        ]
        for mode_trans, mode_rot in itertools.product(modes_trans, modes_rot):
            try:
                propagation = mode_trans | mode_rot
                rotation = (mode_rot != Propagation.NONE) or \
                           (mode_trans == Propagation.SYSTEM_DEFAULT)
                if propagation & (Propagation.SYSTEM_DEFAULT |
                                  Propagation.TRANS_STOKESIAN |
                                  Propagation.ROT_STOKESIAN):
                    pos = x0 + 100. * np.random.random(3)
                else:
                    pos = x0
                system.part.add(pos=pos, v=v0, omega_lab=o0,
                                ext_force=ext_force, ext_torque=ext_torque,
                                rotation=3 * [rotation],
                                mass=0.5, rinertia=[0.5, 0.5, 0.5],
                                propagation=propagation)
            except BaseException:
                continue
        assert len(system.part) > 3
        system.time = 0.
        for i in range(10):
            system.integrator.run(2**i)
            for p in system.part.all():
                pos = np.copy(p.pos)
                vel = np.copy(p.v)
                rot = np.copy(p.omega_lab)
                ref_pos, ref_vel, ref_rot, atol = calc_trajectory(p)
                np.testing.assert_allclose(pos, ref_pos, rtol=0., atol=atol)
                np.testing.assert_allclose(vel, ref_vel, rtol=0., atol=atol)
                np.testing.assert_allclose(rot, ref_rot, rtol=0., atol=atol)


if __name__ == "__main__":
    ut.main()
