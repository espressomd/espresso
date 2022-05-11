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
import unittest as ut
import unittest_decorators as utx
from lb_stats import TestLB

import espressomd
import espressomd.lb


@ut.skipIf(TestLB.n_nodes > 1,
           "LB with N-square only works on 1 MPI rank")
class TestLBCPU(TestLB, ut.TestCase):

    def setUp(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=0)
        self.lb_class = espressomd.lb.LBFluid
        self.params.update({"mom_prec": 1E-9, "mass_prec_per_node": 5E-8})


@utx.skipIfMissingGPU()
class TestLBGPU(TestLB, ut.TestCase):

    def setUp(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=0)
        self.lb_class = espressomd.lb.LBFluidGPU
        self.params.update({"mom_prec": 1E-3, "mass_prec_per_node": 1E-5})


if __name__ == "__main__":
    ut.main()
