#
# Copyright (C) 2013-2018 The ESPResSo project
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
import espressomd

from espressomd.electrostatics import MMM1D


@utx.skipIfMissingFeatures("ELECTROSTATICS")
class ElectrostaticInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        if not self.system.part:
            self.system.periodicity = [0, 0, 1]
            self.system.cell_system.set_n_square()
            self.system.box_l = [10, 10, 10]
            self.system.part.add(id=0, pos=[0, 0, 0])
            self.system.part.add(id=1, pos=[0.1, 0.1, 0.1])
            self.system.part[0].q = 1
            self.system.part[1].q = -1

    if espressomd.has_features("ELECTROSTATICS"):
        test_mmm1d = tests_common.generate_test_for_class(
            system, MMM1D, dict(prefactor=2.0,
                                maxPWerror=0.001,
                                far_switch_radius=3,
                                tune=False))


if __name__ == "__main__":
    ut.main()
