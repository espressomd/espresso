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


class SymplecticEuler(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        self.system.non_bonded_inter.reset()

    def test_integrator_activation(self):
        """Test that symplectic Euler integrator can be activated."""
        # Set symplectic Euler integrator
        self.system.integrator.set_symplectic_euler()

        # Should not raise any exception
        p = self.system.part.add(pos=[0, 0, 0], v=[1, 0, 0], mass=1.0)
        self.system.integrator.run(1)

        # Particle should have moved
        self.assertNotEqual(p.pos[0], 0.0)

    def test_free_particle_propagation(self):
        """Test free particle propagation with symplectic Euler."""
        # Add a free particle
        v_init = [1., 2., -3.5]
        p = self.system.part.add(pos=[0, 0, 0], v=v_init)

        self.system.integrator.set_symplectic_euler()

        dt = self.system.time_step
        expected_pos = np.copy(3. * dt * p.v)

        self.system.integrator.run(3)

        # Check position
        np.testing.assert_allclose(np.copy(p.pos), expected_pos)

        # Check velocity (should remain unchanged for free particle)
        np.testing.assert_allclose(np.copy(p.v), v_init)

    def test_velocity_update_order(self):
        """Test that velocity is updated before position in symplectic Euler."""
        # Add particle with force (use external force for simplicity)
        p = self.system.part.add(pos=[0, 0, 0], v=[0, 0, 0], mass=2.)
        p.ext_force = [1.0, 2.5, -1.]  # Constant external force

        self.system.integrator.set_symplectic_euler()

        # In symplectic Euler:
        # v(n+1) = v(n) + dt * F(n) / m
        # x(n+1) = x(n) + dt * v(n+1)

        dt = self.system.time_step
        expected_v = np.copy(dt * p.ext_force / p.mass)
        expected_x = np.copy(dt * expected_v)  # dt * v(n+1)

        self.system.integrator.run(1)

        # Check that velocity was updated first, then position
        np.testing.assert_allclose(np.copy(p.v), expected_v, atol=1E-10)
        np.testing.assert_allclose(np.copy(p.pos), expected_x, atol=1E-10)


if __name__ == "__main__":
    ut.main()
