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
from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from espressomd.interactions import FeneBond
from espressomd.utils import handle_errors

from tests_common import verify_lj_forces
from numpy import random
from virtual_sites_tracers_common import VirtualSitesTracersCommon


required_features ="VIRTUAL_SITES_INERTIALESS_TRACERS","LB"
@ut.skipIf(not espressomd.has_features(required_features),
           "Test requires VIRTUAL_SITES_INERTIALESS_TRACERS")
class VirtualSitesTracers(ut.TestCase,VirtualSitesTracersCommon):
    if espressomd.has_features(required_features):
        box_height = 10. 
        box_lw=8.
        system = espressomd.System(box_l=(box_lw,box_lw,box_height))
        system.time_step = 0.03
        system.cell_system.skin = 0.1
        lbf = lb.LBFluid(
            agrid=1, dens=1, visc=2, tau= system.time_step, fric = 1)
        system.actors.add(lbf)
        system.thermostat.set_lb(kT=0,act_on_virtual=False)
           
        # Setup boundaries
        walls = [lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=shapes.Wall(normal=[0,0,1], dist = 0.5))
        walls[1].set_params(shape=shapes.Wall(normal=[0,0,-1], dist = -box_height-0.5))
            
        for wall in walls:
           system.lbboundaries.add(wall)
            
        handle_errors("setup") 
    
if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
