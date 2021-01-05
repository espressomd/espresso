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
import mmm1d
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["ELECTROSTATICS", "MMM1D_GPU"])
@utx.skipIfMissingGPU()
class MMM1D_GPU_Test(mmm1d.ElectrostaticInteractionsTests, ut.TestCase):

    def setUp(self):
        self.MMM1D = espressomd.electrostatics.MMM1DGPU
        super().setUp()


if __name__ == "__main__":
    ut.main()
