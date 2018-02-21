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
from espressomd.utils import handle_errors

from tests_common import verify_lj_forces
from numpy import random


class VirtualSitesTracersCommon(object):
    
    def test_aa_method_switching(self):
        # Virtual sites should be disabled by default
        self.assertTrue(isinstance(self.system.virtual_sites, VirtualSitesOff))

        # Switch implementation
        self.system.virtual_sites=VirtualSitesInertialessTracers()
        self.assertTrue(isinstance(self.system.virtual_sites, VirtualSitesInertialessTracers))
        self.assertEqual(self.system.virtual_sites.have_velocity,True)


    def test_advection(self):

        # System setup
        system = self.system
        box_lw=self.box_lw
        box_height = self.box_height

        system.virtual_sites=VirtualSitesInertialessTracers()
        self.lbf.set_params(ext_force=(0.1,0.,0.)) 
        
        
        # Establish steady state flow field
        system.part.add(id=0, pos=(0,5.5,5.5),virtual=1)
        system.integrator.run(500)
            
        # 
        # 
        system.part[0].pos=(0,5.5,5.5)
        
        system.time=0
        ## Perform integration
        
        print("time, actual position, expected position")
        for i in range(10):
            system.integrator.run(100)
            # compute expected position
            X = self.lbf[0,5,5].velocity[0]*system.time
            print( system.time, system.part[0].pos[0],X,system.part[0].pos[0]-X)
            self.assertAlmostEqual(system.part[0].pos[0]/X-1,0,delta=0.005)
        

    def stop_fluid(self):
        system=self.system
        for i in range(int(system.box_l[0])):
            for j in range(int(system.box_l[1])):
                for k in range(int(system.box_l[2])):
                   self.lbf[i,j,k].velocity=(0.,0.,0.)
    
    def compute_angle(self):
        system=self.system
        pos0 = system.part[0].pos
        pos1 = system.part[1].pos
        pos2 = system.part[2].pos
        pos3 = system.part[3].pos    
        
        
        # first normal vector
        n1 = np.cross( (pos1-pos0), (pos2-pos0))
        n2 = np.cross( (pos2-pos0), (pos3-pos0))
    
        
        norm1 = np.linalg.norm(n1)
        norm2 = np.linalg.norm(n2)
        n1 = n1 / norm1
        n2 = n2 / norm2
        
        
        alpha = np.arccos ( np.dot(n1, n2) )
        return alpha
    
    @ut.skipIf(not espressomd.has_features("IMMERSED_BOUNDARY"), "skipped for lack of IMMERSED_BOUNDARY")
    def test_tribend(self):
              
        
        # two triangles with bending interaction
        # move nodes, should relax back

        system=self.system
        system.virtual_sites=VirtualSitesInertialessTracers()
        self.lbf.set_params(ext_force=(0.0,0.,0.)) 
        self.stop_fluid()
        
        system.part.clear()

        ## Add four particles
        system.part.add(id=0, pos=[5, 5, 5], virtual=1)    
        system.part.add(id=1, pos=[5, 5, 6], virtual=1)
        system.part.add(id=2, pos=[5, 6, 6], virtual=1)
        system.part.add(id=3, pos=[5, 6, 5], virtual=1)
        
        ## Add first triel, weak modulus
        from espressomd.interactions import IBM_Triel
        tri1 = IBM_Triel(ind1=0, ind2=1, ind3=2, elasticLaw="Skalak", k1=0.1, k2=0, maxDist = 2.4)
        system.bonded_inter.add(tri1)
        system.part[0].add_bond((tri1, 1, 2))
        
        ## Add second triel
        tri2 = IBM_Triel(ind1=0, ind2=2, ind3=3, elasticLaw="Skalak", k1=10, k2=0, maxDist = 2.4)
        system.bonded_inter.add(tri2)
        system.part[0].add_bond((tri2, 2, 3))
        
        ## Add bending
        from espressomd.interactions import IBM_Tribend
        tribend = IBM_Tribend(ind1=0, ind2=1, ind3=2,ind4=3,kb=1, refShape = "Initial")
        system.bonded_inter.add(tribend)
        system.part[0].add_bond((tribend, 1, 2, 3))
        
        ## output before
        print( "Angle before twisting: " + str(self.compute_angle()))
        
        ## twist
        system.part[1].pos = [5.2, 5, 6]
        
        ## output after
        print( "Angle after twisting: " + str(self.compute_angle()))
        
        ## Perform integrat[ion
        last_angle=self.compute_angle()
        for i in range(8):
            system.integrator.run(500)
            angle=self.compute_angle()
            print( "Angle after relaxation: ",angle)
            self.assertLess(angle,last_angle)
            last_angle=angle

        self.assertLess(angle,0.03)
    
    @ut.skipIf(not espressomd.has_features("IMMERSED_BOUNDARY"),"skipped for lack of IMMERSED_BOUNDARY")
    def test_triel(self):
        system=self.system
        system.virtual_sites=VirtualSitesInertialessTracers()
        self.lbf.set_params(ext_force=(0.0,0.,0.)) 
        self.stop_fluid()
        system.virtual_sites=VirtualSitesInertialessTracers()
        #system.integrator.run(1000)
        
        
        system.part.clear()
        ## Add particles: 0-2 are non-bonded, 3-5 are weakly bonded, 6-8 are strongly bonded
        system.part.add(id=0, pos=[5, 5, 5], virtual=1)    
        system.part.add(id=1, pos=[5, 5, 6], virtual=1)
        system.part.add(id=2, pos=[5, 6, 6], virtual=1)
        
        system.part.add(id=3, pos=[2, 5, 5], virtual=1)    
        system.part.add(id=4, pos=[2, 5, 6], virtual=1)
        system.part.add(id=5, pos=[2, 6, 6], virtual=1)
        
        system.part.add(id=6, pos=[4, 7, 7], virtual=1)    
        system.part.add(id=7, pos=[4, 7, 8], virtual=1)
        system.part.add(id=8, pos=[4, 8, 8], virtual=1)
        
        ## Add triel, weak modulus for 3-5
        from espressomd.interactions import IBM_Triel
        triWeak = IBM_Triel(ind1=3, ind2=4, ind3=5, elasticLaw="Skalak", k1=5, k2=0, maxDist = 2.4)
        system.bonded_inter.add(triWeak)
        system.part[3].add_bond((triWeak, 4, 5))
        
        ## Add triel, strong modulus for 6-8
        triStrong = IBM_Triel(ind1=6, ind2=7, ind3=8, elasticLaw="Skalak", k1=15, k2=0, maxDist = 2.4)
        system.bonded_inter.add(triStrong)
        system.part[6].add_bond((triStrong, 7, 8))
        ## Perform integration
        system.integrator.run(1)
        system.part[3:].pos =system.part[3:].pos +random.random((6,3)) -.5
        
        
        system.integrator.run(15000)
        # get new shapes
        dist1non = np.linalg.norm( np.array( system.part[1].pos - system.part[0].pos ) )
        dist2non = np.linalg.norm( np.array( system.part[2].pos - system.part[0].pos ) )
    
        dist1weak = np.linalg.norm( np.array( system.part[3].pos - system.part[4].pos ) )
        dist2weak = np.linalg.norm( np.array( system.part[3].pos - system.part[5].pos ) )
    
        dist1strong = np.linalg.norm( np.array( system.part[6].pos - system.part[7].pos ) )
        dist2strong = np.linalg.norm( np.array( system.part[6].pos - system.part[8].pos ) )
    
    
        print( "** Distances: non-bonded, weak, strong, expected")
        print( str(dist1non) + "    " + str(dist1weak) + "     " + str(dist1strong) + "    1")
        print( str(dist2non) + "    " + str(dist2weak) + "     " + str(dist2strong) + "    1.414")
        #self.assertAlmostEqual(dist1non,1,delta=0.03)
        self.assertAlmostEqual(dist1weak,1,delta=0.03)
        self.assertAlmostEqual(dist1strong,1,delta=0.03)
        
        #self.assertAlmostEqual(dist2non,np.sqrt(2),delta=0.03)
        self.assertAlmostEqual(dist2weak,np.sqrt(2),delta=0.03)
        self.assertAlmostEqual(dist2strong,np.sqrt(2),delta=0.03)

    
        
