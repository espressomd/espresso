# Copyright (C) 2010-2022 The ESPResSo project
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

import espressomd
import espressomd.lb
import unittest as ut
import unittest_decorators as utx
from engine_lb import SwimmerTest


@utx.skipIfMissingFeatures(["ENGINE", "ROTATIONAL_INERTIA", "MASS"])
@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
class SwimmerTestNSquareCPU(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluid
    tol = 1e-10

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["ENGINE", "ROTATIONAL_INERTIA", "MASS"])
@ut.skipIf(SwimmerTest.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
class SwimmerTestNSquareGPU(SwimmerTest, ut.TestCase):

    lb_class = espressomd.lb.LBFluidGPU
    tol = 1e-5

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


if __name__ == "__main__":
    ut.main()
