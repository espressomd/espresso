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
import espressomd.propagation
import numpy as np
import unittest as ut


class SymplecticEuler(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = 0.01

    def calc_pos(self, p, x0, v0):
        t = self.system.time
        return np.copy(0.5 * p.ext_force / p.mass * t**2 + v0 * t + x0)

    def test_newton_laws(self):
        """
        Tests integration of Newton's equations for a particle with and without
        external forces and with time step changes on the way.

        """
        system = self.system
        system.integrator.set_symplectic_euler()

        # Newton's 1st law with time step change on the way
        p = system.part.add(pos=(0, 0, 0), v=(1, 2, 3))
        system.time_step = 0.02
        system.time = 12.
        total_steps = 100
        np.testing.assert_allclose(np.copy(p.v), (1, 2, 3))
        for i in range(total_steps):
            np.testing.assert_allclose(np.copy(p.pos), np.copy(
                i * system.time_step * p.v), atol=1E-12)
            system.integrator.run(1)

        # Check that the time has passed
        np.testing.assert_allclose(
            system.time, 12. + total_steps * system.time_step)

        # Newton's 2nd law
        if espressomd.has_features("EXTERNAL_FORCES"):
            if espressomd.has_features("MASS"):
                p.mass = 2.3
            max_local_error = list()
            for t in range(4):
                system.time = 0.
                p.pos = (0, 0, 0)
                pos0 = np.copy(p.pos)
                p.v = (0, 0, 0)
                v0 = np.copy(p.v)
                ext_force = np.array([-2., 1.3, 1.])
                p.ext_force = ext_force
                system.time_step = 0.001 * 2**t
                local_error = list()
                for i in range(total_steps):
                    ref_pos = self.calc_pos(p, pos0, v0)
                    np.testing.assert_allclose(
                        np.copy(p.pos), ref_pos, rtol=1E-1, atol=1E-4)
                    local_error.append(max(abs(p.pos - ref_pos)))
                    system.integrator.run(1)
                max_local_error.append(max(local_error))

            # Maximum local error of position within symplectic Euler is
            # approximately proportional to dt^2
            for i in range(3):
                self.assertAlmostEqual(
                    max_local_error[i + 1] / max_local_error[i], 4.0, delta=0.4)


if __name__ == "__main__":
    ut.main()
