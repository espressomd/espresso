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
import espressomd.magnetostatic_extensions


@utx.skipIfMissingFeatures(["DIPOLES"])
class MagnetostaticsInterface(ut.TestCase):
    system = espressomd.System(box_l=[10., 10., 10.])

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.periodicity = [True, True, True]
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
        self.system.actors.clear()
        actor_dawaanr = espressomd.magnetostatics.DipolarDirectSumCpu(
            prefactor=1)
        actor_mdlc = espressomd.magnetostatic_extensions.DLC(
            gap_size=2, maxPWerror=1e-5)
        self.system.actors.add(actor_dawaanr)
        with self.assertRaisesRegex(RuntimeError, 'MDLC cannot extend the currently active magnetostatics solver'):
            self.system.actors.add(actor_mdlc)
        self.system.actors.clear()
        actor = espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu(
            prefactor=1, n_replica=-2)
        with self.assertRaisesRegex(RuntimeError, 'requires n_replica >= 0'):
            self.system.actors.add(actor)
        self.system.actors.clear()
        actor = espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu(
            prefactor=1, n_replica=0)
        with self.assertRaisesRegex(RuntimeError, 'with replica does not support a periodic system with zero replica'):
            self.system.actors.add(actor)
        self.system.actors.clear()
        self.system.periodicity = [True, True, False]
        actor = espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu(
            prefactor=1, n_replica=1)
        solver_mdlc = espressomd.magnetostatic_extensions.DLC(
            gap_size=1, maxPWerror=1e-5)
        self.system.actors.add(actor)
        with self.assertRaisesRegex(RuntimeError, "MDLC requires periodicity 1 1 1"):
            self.system.actors.add(solver_mdlc)
        self.system.actors.remove(solver_mdlc)
        self.system.periodicity = [True, True, True]
        self.system.box_l = [10., 10. + 2e-3, 10.]
        with self.assertRaisesRegex(RuntimeError, "box size in x direction is different from y direction"):
            self.system.actors.add(solver_mdlc)
        self.system.actors.clear()
        self.system.box_l = [10., 10., 10.]


if __name__ == "__main__":
    ut.main()
