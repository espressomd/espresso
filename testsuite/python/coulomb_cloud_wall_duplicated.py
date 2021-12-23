#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd.cuda_init
import espressomd.electrostatics
import tests_common


@utx.skipIfMissingFeatures("ELECTROSTATICS")
class CoulombCloudWall(ut.TestCase):

    """This compares p3m, p3m_gpu electrostatic forces and energy against
    stored data."""
    system = espressomd.System(box_l=[10., 10., 20.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    data = np.genfromtxt(
        tests_common.abspath("data/coulomb_cloud_wall_duplicated_system.data"))

    tolerance = 1E-3

    # Reference energy from p3m in the tcl test case
    reference_energy = 2. * 148.94229549

    def setUp(self):
        # Add particles to system and store reference forces in hash
        # Input format: id pos q f
        self.system.part.add(pos=self.data[:, 1:4], q=self.data[:, 4])
        self.forces = self.data[:, 5:8]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def compare(self, method_name, energy=True):
        # Compare forces and energy now in the system to stored ones

        # Force
        force_diff = np.linalg.norm(
            self.system.part.all().f - self.forces, axis=1)
        self.assertLess(
            np.mean(force_diff), self.tolerance,
            msg="Absolute force difference too large for method " + method_name)

        # Energy
        if energy:
            self.assertAlmostEqual(
                self.system.analysis.energy()["total"], self.reference_energy,
                delta=self.tolerance,
                msg="Absolute energy difference too large for " + method_name)

    # Tests for individual methods

    @utx.skipIfMissingFeatures("P3M")
    def test_p3m(self):
        self.system.actors.add(
            espressomd.electrostatics.P3M(
                prefactor=1, r_cut=1.001, accuracy=1e-3,
                mesh=[64, 64, 128], cao=7, alpha=2.70746, tune=False))
        self.system.integrator.run(0)
        self.compare("p3m", energy=True)

    @utx.skipIfMissingGPU()
    def test_p3m_gpu(self):
        self.system.actors.add(
            espressomd.electrostatics.P3MGPU(
                prefactor=1,
                r_cut=1.001,
                accuracy=1e-3,
                mesh=[64, 64, 128],
                cao=7,
                alpha=2.70746,
                tune=False))
        self.system.integrator.run(0)
        self.compare("p3m_gpu", energy=False)

    def test_zz_deactivation(self):
        # The energy is 0 if no method is active
        self.assertEqual(self.system.analysis.energy()["total"], 0.0)


if __name__ == "__main__":
    ut.main()
