
#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function
import espressomd
import numpy as np
import unittest as ut


class Integrate(ut.TestCase):

    """
    Tests integration of Newton's equations for a particle with and without
    external forces and with time step changes on the way.

    """

    def test(self):
        system = espressomd.System(box_l=[10.0, 10.0, 10.0])
        system.cell_system.skin = 0

        # Newton's 1st law with time step change on the way
        p = system.part.add(pos=(0, 0, 0), v=(1, 2, 3))
        system.time_step = 0.01
        system.time = 12.
        np.testing.assert_allclose(np.copy(p.v), (1, 2, 3))
        for i in range(10):
            np.testing.assert_allclose(np.copy(p.pos), np.copy(
                i * system.time_step * p.v), atol=1E-12)
            system.integrator.run(1)

        # Check that the time has passed
        np.testing.assert_allclose(system.time, 12. + 10 * system.time_step)

        v = p.v
        pos1 = p.pos
        system.time_step = 0.02
        np.testing.assert_allclose(np.copy(v), np.copy(p.v), atol=1E-12)
        np.testing.assert_allclose(np.copy(pos1), np.copy(p.pos), atol=1E-12)
        for i in range(10):
            np.testing.assert_allclose(np.copy(p.pos), np.copy(
                pos1 + i * system.time_step * p.v), atol=1E-12)
            system.integrator.run(1)

        # Newton's 2nd law
        if espressomd.has_features("EXTERNAL_FORCES"):
            if espressomd.has_features("MASS"):
                p.mass = 2.3
            p.pos = (0, 0, 0)
            ext_force = np.array((-2, 1.3, 1))
            p.ext_force = ext_force
            system.time_step = 0.03
            for i in range(10):
                np.testing.assert_allclose(np.copy(p.pos), np.copy(
                    0.5 * ext_force / p.mass * (i * system.time_step)**2 + v * i * system.time_step), atol=1E-12)
                system.integrator.run(1)


if __name__ == '__main__':
    ut.main()
