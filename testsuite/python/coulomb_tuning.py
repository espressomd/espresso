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
import espressomd.cuda_init
import espressomd.electrostatics
import tests_common


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class CoulombCloudWallTune(ut.TestCase):

    """This compares p3m, p3m_gpu electrostatic forces against stored data."""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    tolerance = 1E-3

    def setUp(self):
        self.system.box_l = (10, 10, 10)
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4

        data = np.load(tests_common.abspath("data/coulomb_tuning_system.npz"))
        self.forces = data['forces']
        self.system.part.add(pos=data['pos'], q=data['charges'])

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()

    def compare(self, method_name):
        # Compare forces now in the system to stored ones
        difference = np.linalg.norm(
            self.system.part.all().f - self.forces, axis=1)
        self.assertLessEqual(
            np.mean(difference), self.tolerance,
            "Absolute force difference too large for method " + method_name)

    # Tests for individual methods
    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m(self):
        # We have to add some tolerance here, because the reference
        # system is not homogeneous
        self.system.actors.add(
            espressomd.electrostatics.P3M(prefactor=1., accuracy=5e-4,
                                          tune=True))
        self.system.integrator.run(0)
        self.compare("p3m")

    @utx.skipIfMissingGPU()
    def test_p3m_gpu(self):
        # We have to add some tolerance here, because the reference
        # system is not homogeneous
        self.system.actors.add(
            espressomd.electrostatics.P3MGPU(prefactor=1., accuracy=5e-4,
                                             tune=True))
        self.system.integrator.run(0)
        self.compare("p3m_gpu")


if __name__ == "__main__":
    ut.main()
