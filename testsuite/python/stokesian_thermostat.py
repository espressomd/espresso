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
import numpy as np


@utx.skipIfMissingFeatures("STOKESIAN_DYNAMICS")
class StokesianThermostat(ut.TestCase):

    """Test Stokesian thermostat"""
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
        with self.assertRaises(RuntimeError):
            system.thermostat.set_stokesian(kT=-1)
        with self.assertRaises(ValueError):
            system.thermostat.set_stokesian(kT=1, seed=-1)

        system.thermostat.set_stokesian(kT=kT, seed=41)
        self.system.integrator.set_stokesian_dynamics(
            viscosity=viscosity, radii={0: radius})

        pos2force = np.sqrt(24 * kT * system.time_step)

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

        # Same seed should not give the same force with different counter state
        # force1: brownian.rng_counter() = 0, brownian.rng_seed() = 41
        # force2: brownian.rng_counter() = 1, brownian.rng_seed() = 41
        p = reset_particle()
        system.thermostat.set_stokesian(kT=kT, seed=41)
        system.integrator.run(1)
        force2 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force2, force1)))

        # Seed offset should not give the same force with a lag
        # force3: brownian.rng_counter() = 2, brownian.rng_seed() = 42
        # force4: brownian.rng_counter() = 3, brownian.rng_seed() = 41
        p = reset_particle()
        system.thermostat.set_stokesian(kT=kT, seed=42)
        system.integrator.run(1)
        force3 = np.copy(p.pos) / pos2force
        p = reset_particle()
        system.thermostat.set_stokesian(kT=kT, seed=41)
        system.integrator.run(1)
        force4 = np.copy(p.pos) / pos2force
        self.assertTrue(np.all(np.not_equal(force3, force4)))

    def test_integrator_exceptions(self):
        # invalid parameters should throw exceptions
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=1.0, radii={0: -1})
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=-1, radii={0: 1.0})

        # invalid PBC should throw exceptions
        self.system.integrator.set_vv()
        self.system.periodicity = [False, False, True]
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_stokesian_dynamics(
                viscosity=1.0, radii={0: 1.0})

        self.system.periodicity = [False, False, False]
        self.system.integrator.set_stokesian_dynamics(
            viscosity=1.0, radii={0: 1.0})

        with self.assertRaises(Exception):
            self.system.periodicity = [False, True, False]


if __name__ == "__main__":
    ut.main()
