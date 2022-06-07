#
# Copyright (C) 2017-2019 The ESPResSo project
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
        self.system.actors.clear()
        self.system.part.clear()

    def compare(self, actor):
        self.system.actors.add(actor)
        self.system.integrator.run(0)
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f), self.ref_forces, atol=2e-3)

    def test_p3m_cpu(self):
        actor = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=5e-4, tune=True)
        self.compare(actor)

    @utx.skipIfMissingGPU()
    def test_p3m_gpu(self):
        actor = espressomd.electrostatics.P3MGPU(
            prefactor=1., accuracy=5e-4, tune=True)
        self.compare(actor)


if __name__ == "__main__":
    ut.main()
