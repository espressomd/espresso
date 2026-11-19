#
# Copyright (C) 2017-2026 The ESPResSo project
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
import numpy as np
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.electrostatics
import tests_common


@utx.skipIfMissingFeatures(["P3M"])
class CoulombCloudWallTune(ut.TestCase):

    """This compares P3M electrostatic forces against reference values."""
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        data = np.load(tests_common.data_path("coulomb_tuning_system.npz"))
        self.ref_forces = data['forces']
        self.system.part.add(pos=data['pos'], q=data['charges'])

    def tearDown(self):
        self.system.electrostatics.clear()
        self.system.part.clear()

    def compare(self, actor):
        self.system.electrostatics.solver = actor
        self.system.integrator.run(0)
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f), self.ref_forces, atol=2e-3)

    def test_p3m_cpu(self):
        actor = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=5e-4, tune=True)
        self.compare(actor)

    def test_p3m_cpu_tune_static(self):
        # default: timings are force calculations on fixed positions
        pos_before = np.copy(self.system.part.all().pos)
        actor = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=5e-4, tune=True, timings=2)
        self.assertFalse(actor.tune_with_integration)
        self.compare(actor)
        np.testing.assert_allclose(
            np.copy(self.system.part.all().pos), pos_before, rtol=0., atol=0.)

    def test_p3m_cpu_tune_with_integration(self):
        # opt-in: timings are full integration steps, the particles move
        pos_before = np.copy(self.system.part.all().pos)
        actor = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=5e-4, tune=True, timings=2,
            tune_with_integration=True)
        self.assertTrue(actor.tune_with_integration)
        self.system.electrostatics.solver = actor
        self.assertTrue(actor.is_tuned)
        displacement = np.linalg.norm(
            np.copy(self.system.part.all().pos) - pos_before, axis=1)
        self.assertGreater(np.max(displacement), 0.)
        # the tuned parameters still reach the requested accuracy: restore
        # the reference configuration and check the forces
        self.system.part.all().pos = pos_before
        self.system.integrator.run(0)
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f), self.ref_forces, atol=2e-3)

    def test_p3m_cpu_tune_skin(self):
        # opt-in: the skin is tuned for the fastest parameter sets
        pos_before = np.copy(self.system.part.all().pos)
        skin_before = self.system.cell_system.skin
        actor = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=5e-4, tune=True, timings=2,
            tune_with_integration=True, tune_skin=True)
        self.assertTrue(actor.tune_skin)
        self.system.electrostatics.solver = actor
        self.assertTrue(actor.is_tuned)
        skin = self.system.cell_system.skin
        self.assertNotAlmostEqual(skin, skin_before)
        self.assertGreater(skin, 0.)
        # the tuned skin must fit in the cell system together with r_cut
        self.assertLess(actor.r_cut + skin, self.system.box_l[0] / 2.)
        # the tuned parameters still reach the requested accuracy
        self.system.part.all().pos = pos_before
        self.system.integrator.run(0)
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f), self.ref_forces, atol=2e-3)
        self.system.cell_system.skin = skin_before

    @utx.skipIfMissingGPU()
    def test_p3m_gpu(self):
        actor = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=5e-4, tune=True, gpu=True)
        self.compare(actor)


if __name__ == "__main__":
    ut.main()
