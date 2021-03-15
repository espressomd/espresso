#
# Copyright (C) 2013-2020 The ESPResSo project
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
import tests_common

import espressomd.magnetostatics


@utx.skipIfMissingFeatures(["DIPOLES"])
class MagnetostaticsInterface(ut.TestCase):
    system = espressomd.System(box_l=[10., 10., 10.])

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.part.add(pos=(0.0, 0.0, 0.0), dip=(1.3, 2.1, -6))
        self.system.part.add(pos=(0.1, 0.1, 0.1), dip=(7.3, 6.1, -4))

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    if espressomd.has_features("DIPOLES"):
        test_dds_cpu = tests_common.generate_test_for_class(
            system, espressomd.magnetostatics.DipolarDirectSumCpu,
            dict(prefactor=3.4))

    if espressomd.has_features(
            "DIPOLAR_DIRECT_SUM") and espressomd.gpu_available():
        test_dds_gpu = tests_common.generate_test_for_class(
            system, espressomd.magnetostatics.DipolarDirectSumGpu,
            dict(prefactor=3.4))

    if espressomd.has_features("DIPOLES"):
        test_dds_replica = tests_common.generate_test_for_class(
            system, espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu,
            dict(prefactor=3.4, n_replica=2))

    def test_exceptions(self):
        actor = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=-1)
        with self.assertRaises(ValueError):
            self.system.actors.add(actor)
        actor = self.system.actors[0]
        actor.set_params(prefactor=1)
        with self.assertRaises(ValueError):
            actor.set_params(prefactor=-1)
        with self.assertRaises(ValueError):
            actor.set_magnetostatics_prefactor()


if __name__ == "__main__":
    ut.main()
