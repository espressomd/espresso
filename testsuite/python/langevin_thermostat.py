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


class LangevinThermostat(ut.TestCase):

    """Test Langevin Dynamics"""
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
            p = system.part.add(pos=[0, 0, 0])
            return p

        system = self.system
        system.time_step = 0.01

        kT = 1.1
        gamma = 3.5

        # No seed should throw exception
        with self.assertRaises(ValueError):
            system.thermostat.set_langevin(kT=kT, gamma=gamma)

        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=41)

        # run(0) does not increase the philox counter and should give the same
        # force
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force0 = np.copy(p.f)
        system.integrator.run(0, recalc_forces=True)
        force1 = np.copy(p.f)
        np.testing.assert_almost_equal(force0, force1)

        # run(1) should give a different force
        p = reset_particle()
        system.integrator.run(1)
        force2 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force1, force2)))

        # Different seed should give a different force with same counter state
        # force2: langevin.rng_counter() = 1, langevin.rng_seed() = 41
        # force3: langevin.rng_counter() = 1, langevin.rng_seed() = 42
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force2 = np.copy(p.f)
        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
        system.integrator.run(0, recalc_forces=True)
        force3 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force2, force3)))

        # Same seed should not give the same force with different counter state
        p = reset_particle()
        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
        system.integrator.run(1)
        force4 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force3, force4)))

        # Seed offset should not give the same force with a lag
        # force4: langevin.rng_counter() = 2, langevin.rng_seed() = 42
        # force5: langevin.rng_counter() = 3, langevin.rng_seed() = 41
        reset_particle()
        system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=41)
        system.integrator.run(1)
        p = reset_particle()
        system.integrator.run(0, recalc_forces=True)
        force5 = np.copy(p.f)
        self.assertTrue(np.all(np.not_equal(force4, force5)))

    def test_02__friction_trans(self):
        """Tests the translational friction-only part of the thermostat."""

        system = self.system
        # Translation
        gamma_t_i = 2
        gamma_t_a = np.array((0.5, 2, 1.5))
        v0 = np.array((5., 5., 5.))

        system.time_step = 0.0005
        system.part.add(pos=(0, 0, 0), v=v0)
        if espressomd.has_features("MASS"):
            system.part[0].mass = 3
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            system.thermostat.set_langevin(kT=0, gamma=gamma_t_a, seed=41)
        else:
            system.thermostat.set_langevin(kT=0, gamma=gamma_t_i, seed=41)

        system.time = 0
        for _ in range(100):
            system.integrator.run(10)
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                np.testing.assert_allclose(
                    np.copy(system.part[0].v),
                    v0 * np.exp(-gamma_t_a /
                                system.part[0].mass * system.time),
                    atol=4E-4)
            else:
                np.testing.assert_allclose(
                    np.copy(system.part[0].v),
                    v0 * np.exp(-gamma_t_i /
                                system.part[0].mass * system.time),
                    atol=45E-4)

    @utx.skipIfMissingFeatures("ROTATION")
    def test_03__friction_rot(self):
        """Tests the rotational friction-only part of the thermostat."""

        system = self.system
        # Translation
        gamma_t_i = 2
        gamma_t_a = [0.5, 2, 1.5]
        gamma_r_i = 3
        gamma_r_a = np.array((1.5, 0.7, 1.2))
        o0 = np.array((5., 5., 5.))

        system.time_step = 0.0001
        system.part.add(pos=(0, 0, 0), omega_body=o0, rotation=(1, 1, 1))
        if espressomd.has_features("ROTATIONAL_INERTIA"):
            system.part[0].rinertia = [2, 2, 2]
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            system.thermostat.set_langevin(
                kT=0, gamma=gamma_t_a, gamma_rotation=gamma_r_a, seed=41)
        else:
            system.thermostat.set_langevin(
                kT=0, gamma=gamma_t_i, gamma_rotation=gamma_r_i, seed=41)

        system.time = 0
        if espressomd.has_features("ROTATIONAL_INERTIA"):
            rinertia = np.copy(system.part[0].rinertia)
        else:
            rinertia = np.array((1, 1, 1))
        for _ in range(100):
            system.integrator.run(10)
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                np.testing.assert_allclose(
                    np.copy(system.part[0].omega_body),
                    o0 * np.exp(-gamma_r_a / rinertia * system.time), atol=5E-4)
            else:
                np.testing.assert_allclose(
                    np.copy(system.part[0].omega_body),
                    o0 * np.exp(-gamma_r_i / rinertia * system.time), atol=5E-4)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES")
    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01

        virtual = system.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        system.thermostat.set_langevin(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=False, seed=41)

        system.integrator.run(0)

        np.testing.assert_almost_equal(np.copy(virtual.f), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])

        system.thermostat.set_langevin(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=True, seed=41)
        system.integrator.run(0)

        np.testing.assert_almost_equal(np.copy(virtual.f), [-1, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1, 0, 0])


if __name__ == "__main__":
    ut.main()
