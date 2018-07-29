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
from espressomd.virtual_sites import VirtualSitesInertialessTracers, VirtualSitesOff

from tests_common import verify_lj_forces
from numpy import random


#@ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_INERTIALESS_TRACERS","LB"),
#           "Test requires VIRTUAL_SITES_INERTIALESS_TRACERS")
class VirtualSitesTracers(ut.TestCase):
    s = espressomd.System(box_l= [1,1,1])
    s.seed = range(s.cell_system.get_state()["n_nodes"])
    
    def test_advection(self):

        # System setup
        system = self.s
        system.virtual_sites=VirtualSitesInertialessTracers()
        system.time_step = 0.02
        system.cell_system.skin = 0.1
        
        box_height = 16. 
        box_lw=16
        system.box_l = box_lw,box_lw,box_height
        
        lbf = lb.LBFluidGPU(agrid=1, dens=1, visc=2, tau= system.time_step, fric = 1)
        system.actors.add(lbf)
        
        system.thermostat.set_lb(kT=0)
        
        # Setup boundaries
        walls = [lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=shapes.Wall(normal=[0,0,1], dist = 0.5))
        walls[1].set_params(shape=shapes.Wall(normal=[0,0,-1], dist = -box_height-0.5))
        
        for wall in walls:
            system.lbboundaries.add(wall)
        
        # Establish steady state flow field
        system.part.add(id=0, pos=(0,5.5,5.5),virtual=0,ext_force=(10,0,0))
        for i in range(100):
            last_t=system.time
            last_x=system.part[0].pos
            system.integrator.run(500)
            print(system.part[0].v,(system.part[0].pos-last_x)/(system.time-last_t))
if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
