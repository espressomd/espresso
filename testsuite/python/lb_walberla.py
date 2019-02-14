
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
from espressomd.lb_walberla import LbWalberla
import numpy as np

@ut.skipIf(not espressomd.has_features("LB_WALBERLA"),"Skipping for LACK of LB_WALBERLA")
class LbWalberlaTest(ut.TestCase):
    visc = 2.
    agrid = .5
    vel = [1.3,2.2,3.1]

    def test(self):
        system = espressomd.System(box_l=[10] * 3)
        lb = LbWalberla(visc=self.visc, agrid=self.agrid)
        self.assertEqual(lb.visc, self.visc)
        self.assertEqual(lb.agrid, self.agrid)
        lb[0,0,0].velocity = self.vel
        self.assertEqual(lb[0,0,0].velocity, self.vel)
    
if __name__ == "__main__":
    # print("Features: ", espressomd.features())
    ut.main()
