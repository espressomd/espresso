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
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["P3M"])
class P3MTuneSkin(ut.TestCase):
    """Tune the P3M parameters together with the Verlet skin."""
    system = espressomd.System(box_l=[10., 10., 10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    # random charges are not relaxed: cap the forces during the integration
    system.force_cap = 5.

    def setUp(self):
        rng = np.random.default_rng(seed=42)
        n_part = 60
        self.system.part.add(pos=rng.random((n_part, 3)) * 10.,
                             q=np.tile([1., -1.], n_part // 2))

    def tearDown(self):
        self.system.electrostatics.clear()
        self.system.part.clear()

    def test_tune_skin(self):
        pos_before = np.copy(self.system.part.all().pos)
        skin_before = self.system.cell_system.skin
        solver = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-2, tune=True, timings=2,
            tune_with_integration=True, tune_skin=True)
        self.assertTrue(solver.tune_skin)
        self.system.electrostatics.solver = solver
        self.assertTrue(solver.is_tuned)
        skin = self.system.cell_system.skin
        self.assertNotAlmostEqual(skin, skin_before)
        self.assertGreater(skin, 0.)
        # the tuned skin must fit in the cell system together with r_cut
        self.assertLess(solver.r_cut + skin, self.system.box_l[0] / 2.)
        # the tuned system can be integrated
        self.system.part.all().pos = pos_before
        self.system.integrator.run(2)

    def test_tune_without_skin(self):
        skin_before = self.system.cell_system.skin
        solver = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-2, tune=True, timings=2,
            tune_with_integration=True)
        self.assertFalse(solver.tune_skin)
        self.system.electrostatics.solver = solver
        self.assertEqual(self.system.cell_system.skin, skin_before)


if __name__ == "__main__":
    ut.main()
