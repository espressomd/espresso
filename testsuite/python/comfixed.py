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

import unittest as ut
import espressomd
import espressomd.analyze
import espressomd.lb
import numpy as np


class ComFixed(ut.TestCase):

    np.random.seed(seed=42)

    def com(self, parts):
        return np.average(parts.pos, axis=0, weights=parts.mass)

    def test(self):
        dt = 0.01
        skin = 0.4

        system = espressomd.System(box_l=[10., 10., 10.])
        system.time_step = dt
        system.cell_system.skin = skin

        system.thermostat.set_langevin(kT=1., gamma=0.01, seed=41)

        for i in range(100):
            r = [0.5, 1., 1.] * system.box_l * np.random.random(3)
            v = 3 * [0.]
            # Make sure that id and type gaps work correctly
            system.part.add(id=2 * i, pos=r, v=v, type=2 * (i % 2))
        partcls = system.part.all()

        if espressomd.has_features(["MASS"]):
            # Avoid masses too small for the time step
            partcls.mass = 2. * (0.1 + np.random.random(100))

        com_0 = self.com(partcls)

        system.comfixed.types = [0, 2]

        # Interface check
        self.assertEqual(system.comfixed.types, [2, 0])

        for i in range(10):
            com_i = self.com(partcls)

            for j in range(3):
                self.assertAlmostEqual(com_0[j], com_i[j], places=10)

            system.integrator.run(100)


if __name__ == "__main__":
    ut.main()
