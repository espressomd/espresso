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
import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np


class BrownianThermostat(ut.TestCase):

    """Test Brownian Dynamics"""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.set_domain_decomposition(use_verlet_lists=True)
    system.cell_system.skin = 0
    system.periodicity = [0, 0, 0]

    def setUp(self):
        np.random.seed(42)

    def tearDown(self):
        self.system.time_step = 1e-12
        self.system.cell_system.skin = 0.0
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def test_01__rng(self):
        """Test for RNG consistency."""
        def reset_particle():
            self.system.part.clear()
            p = self.system.part.add(pos=[0, 0, 0])
            if espressomd.has_features("ROTATION"):
                p.rotation = [1, 1, 1]
                # Make sure rinertia does not change diff coeff
                if espressomd.has_features("ROTATIONAL_INERTIA"):
                    p.rinertia = [0.4, 0.4, 0.4]
            return p

        system = self.system
        system.time_step = 0.01

        kT = 1.1
        gamma = 3.5
        pos2force = np.sqrt(2 * kT / gamma * system.time_step)
        omega2torque = np.sqrt(2 * kT / gamma * system.time_step)

        # No seed should throw exception
        with self.assertRaises(ValueError):
            system.thermostat.set_brownian(kT=kT, gamma=gamma)

        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        system.integrator.set_brownian_dynamics()

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

    @utx.skipIfMissingFeatures("VIRTUAL_SITES")
    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01

        virtual = system.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=False, seed=41)
        system.integrator.set_brownian_dynamics()

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [1, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])

        virtual.pos = physical.pos = [0, 0, 0]
        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=True, seed=41)
        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])


if __name__ == "__main__":
    ut.main()
