#
# Copyright (C) 2010-2026 The ESPResSo project
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

    system = espressomd.System(box_l=[10., 10., 10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.system.thermostat.set_langevin(kT=1., gamma=0.01, seed=41)

    def tearDown(self):
        self.system.comfixed.types = []
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def com(self, parts):
        return np.average(parts.pos, axis=0, weights=parts.mass)

    def check_com_conserved(self, fixed_types, particle_types):
        """Set up a system whose particles carry the given ``particle_types``,
        fix the centre of mass of ``fixed_types`` and assert that the COM of
        the fixed particles does not drift during integration."""
        system = self.system

        for i in range(100):
            r = [0.5, 1., 1.] * system.box_l * np.random.random(3)
            v = 3 * [0.]
            # Make sure that id and type gaps work correctly
            system.part.add(id=2 * i, pos=r, v=v,
                            type=particle_types[i % len(particle_types)])

        if espressomd.has_features(["MASS"]):
            # Avoid masses too small for the time step
            system.part.all().mass = 2. * (0.1 + np.random.random(100))

        # COM of the subset of particles whose type is fixed.
        distinct_fixed = set(fixed_types)
        fixed_parts = system.part.select(lambda p: p.type in distinct_fixed)
        com_0 = self.com(fixed_parts)

        system.comfixed.types = fixed_types

        for _ in range(10):
            com_i = self.com(fixed_parts)
            for j in range(3):
                self.assertAlmostEqual(com_0[j], com_i[j], places=10)
            system.integrator.run(100)

    def test(self):
        self.check_com_conserved(fixed_types=[0, 2], particle_types=[0, 2])

        # Interface check
        self.system.comfixed.types = [0, 2]
        self.assertEqual(self.system.comfixed.types, [2, 0])

    def test_duplicate_types(self):
        """Regression test for bug #27: a duplicated entry in the ComFixed
        types list used to leave the internal type->index map with a stored
        index >= map size, producing an out-of-bounds vector access in
        ComFixed::apply() during integration. The duplicate must be handled
        safely (deduplicated): fixing the COM of types ``[1, 1]`` must behave
        exactly like fixing the COM of type ``1`` and conserve it."""
        self.check_com_conserved(fixed_types=[1, 1], particle_types=[1, 3])


if __name__ == "__main__":
    ut.main()
