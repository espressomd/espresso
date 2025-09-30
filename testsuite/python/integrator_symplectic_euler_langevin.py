#
# Copyright (C) 2013-2025 The ESPResSo project
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
import espressomd
import numpy as np
import unittest as ut
import unittest_decorators as utx


class SymplecticEulerLangevin(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    @utx.skipIfMissingFeatures(["MASS"])
    def test_langevin_thermostat_compatibility(self):
        """Test that symplectic Euler works with Langevin thermostat."""
        # Add test particles
        p1 = self.system.part.add(pos=[0, 0, 0], v=[1, 0, 0], mass=1.)
        p2 = self.system.part.add(pos=[2, 0, 0], v=[0, 1, 0], mass=2.)

        # Set Langevin thermostat
        self.system.thermostat.set_langevin(kT=1.0, gamma=0.5, seed=42)

        # Set symplectic Euler integrator (should work)
        self.system.integrator.set_symplectic_euler()

        # Run integration (should not raise exception)
        initial_pos = np.copy([p1.pos, p2.pos])
        self.system.integrator.run(10)
        final_pos = np.copy([p1.pos, p2.pos])

        # Particles should have moved
        self.assertFalse(np.allclose(initial_pos, final_pos))

    def test_no_thermostat_compatibility(self):
        """Test that symplectic Euler works without thermostat."""
        self.system.part.add(pos=[0, 0, 0], v=[1, 0, 0])

        # No thermostat
        self.system.thermostat.turn_off()

        # Set symplectic Euler integrator (should work)
        self.system.integrator.set_symplectic_euler()

        # Run integration (should not raise exception)
        self.system.integrator.run(10)

        # Should work fine
        self.assertTrue(True)

    def test_friction_effects(self):
        """Test that friction effects are applied."""
        # High velocity particle should slow down due to friction
        p = self.system.part.add(pos=[0, 0, 0], v=[10, 0, 0])

        # Set thermostat with high friction, low temperature
        self.system.thermostat.set_langevin(kT=0.1, gamma=5.0, seed=42)
        self.system.integrator.set_symplectic_euler()

        initial_v = np.linalg.norm(p.v)
        self.system.integrator.run(100)
        final_v = np.linalg.norm(p.v)

        # Velocity should decrease due to friction
        self.assertLess(final_v, initial_v)

    def test_stochastic_force(self):
        """Test that stochastic forces are applied."""
        # Particle at rest should start moving due to random forces
        p = self.system.part.add(pos=[0, 0, 0], v=[0, 0, 0])

        # Set thermostat with temperature but no initial velocity
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_symplectic_euler()

        initial_pos = np.copy(p.pos)
        self.system.integrator.run(100)
        final_pos = np.copy(p.pos)

        # Particle should move due to random forces
        displacement = np.linalg.norm(final_pos - initial_pos)
        self.assertGreater(displacement, 0.1)


if __name__ == "__main__":
    ut.main()
