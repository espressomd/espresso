#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
import espressomd
from espressomd.magnetostatics import *
from tests_common import *


class MagnetostaticsInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    def setUp(self):
        self.system.box_l = 10, 10, 10
        if not self.system.part.exists(0):
            self.system.part.add(id=0, pos=(0.1, 0.1, 0.1), dip=(1.3, 2.1, -6))
        if not self.system.part.exists(1):
            self.system.part.add(id=1, pos=(0, 0, 0), dip=(7.3, 6.1, -4))

    if espressomd.has_features(["DP3M"]):
        test_DP3M = generate_test_for_class(
            system, DipolarP3M, dict(
                prefactor=1.0, epsilon=0.0, inter=1000, mesh_off=[
                    0.5, 0.5, 0.5], r_cut=2.4, mesh=[
                    8, 8, 8], cao=1, alpha=12, accuracy=0.01, tune=False))

    if espressomd.has_features(["DIPOLAR_DIRECT_SUM"]):
        test_DdsCpu = generate_test_for_class(
            system, DipolarDirectSumCpu, dict(prefactor=3.4))
        test_DdsRCpu = generate_test_for_class(
            system, DipolarDirectSumWithReplicaCpu, dict(
                prefactor=3.4, n_replica=2))


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
