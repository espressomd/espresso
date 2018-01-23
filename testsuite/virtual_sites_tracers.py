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


@ut.skipIf(not espressomd.has_features("VIRTUAL_SITES" ),
           "Test requires VIRTUAL_SITES")
class VirtualSitesTracers(ut.TestCase):
    s = espressomd.System()
    s.seed = range(s.cell_system.get_state()["n_nodes"])
    
    def test_aa_method_switching(self):
        # Virtual sites should be disabled by default
        self.assertTrue(isinstance(self.s.virtual_sites, VirtualSitesOff))

        # Switch implementation
        self.s.virtual_sites=VirtualSitesInertialessTracers()
        self.assertTrue(isinstance(self.s.virtual_sites, VirtualSitesInertialessTracers))
        self.assertEqual(self.s.virtual_sites.have_velocity,True)


    def test_advection(self):

        # System setup
        system = self.s
        system.virtual_sites=VirtualSitesInertialessTracers()
        system.time_step = 0.05
        system.cell_system.skin = 0.1
        
        box_l = 9. 
        system.box_l = [box_l] * 3
        
        lbf = lb.LBFluid(agrid=1, dens=1, visc=10, tau= system.time_step, ext_force=[0.1, 0, 0], fric = 1)
        system.actors.add(lbf)
        
        system.thermostat.set_lb(kT=0)
        
        # Setup boundaries
        walls = [lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=shapes.Wall(normal=[0,0,1], dist = 0.5))
        walls[1].set_params(shape=shapes.Wall(normal=[0,0,-1], dist = -box_l-0.5))
        
        for wall in walls:
            system.lbboundaries.add(wall)
        
        system.part.add(id=0, pos=[0,box_l/2,box_l/2], virtual=1)    
        # Establish steady state flow field
        system.integrator.run(100)
            
        system.part[0].pos=[0,box_l/2,box_l/2]
        
        system.time=0
        ## Perform integration
        
        print("time, actual position, expected position")
        for i in range(10):
            system.integrator.run(100)
            # compute expected position
            X = lbf[0,int(box_l/2),int(box_l/2)].velocity[0] * system.time
            print( system.time, system.part[0].pos[0],X)
            self.assertAlmostEqual(system.part[0].pos[0],X,delta=0.02)

        
        
            
if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
