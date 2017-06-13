
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
from espressomd.interactions import HarmonicBond,Angle_Harmonic
from espressomd.collision_detection import CollisionMode
import numpy as np

@ut.skipIf(not espressomd.has_features("COLLISION_DETECTION"),"Required features not compiled in")
class CollisionDetection(ut.TestCase):
    """Tests interface and functionality of the collision detection / dynamic binding"""

    s = espressomd.System()

    H = HarmonicBond(k=1,r_0=0.1)
    H2 = HarmonicBond(k=1,r_0=0)
    s.bonded_inter.add(H)
    s.bonded_inter.add(H2)
    s.time_step=0.01
    s.cell_system.skin=0
    s.min_global_cut=0.2

    def test_00_interface_and_defaults(self):
        # Is it off by default
        self.assertEquals(self.s.collision_detection.mode,CollisionMode.off)
        # Make sure params cannot be set individually
        with self.assertRaises(Exception):
            self.s.collision_detection.mode=CollisionMode.bind_centers
        
        # Collision modes should be instances of CollisionMode
        with self.assertRaises(Exception):
            self.s.collision_detection.set_params(mode=0)
        # That should work
        self.s.collision_detection.set_params(mode=CollisionMode.off)
        self.assertEquals(self.s.collision_detection.mode,CollisionMode.off)

    def test_bind_centers(self):
        # Check that it leaves particles alone, wehn off
        self.s.collision_detection.set_params(mode=CollisionMode.off)
        
        self.s.part.clear()
        self.s.part.add(pos=(0,0,0),id=0)
        self.s.part.add(pos=(0.1,0,0),id=1)
        self.s.part.add(pos=(0.1,0.3,0),id=2)
        self.s.integrator.run(0)
        self.assertEqual(self.s.part[0].bonds,())
        self.assertEqual(self.s.part[1].bonds,())
        self.assertEqual(self.s.part[2].bonds,())

        # Check that it cannot be activated 
        self.s.collision_detection.set_params(mode=CollisionMode.bind_centers,distance=0.11,bond_centers=self.H)
        self.s.integrator.run(1,recalc_forces=True)
        bond0=((self.s.bonded_inter[0],1),)
        bond1=((self.s.bonded_inter[0],0),)
        self.assertTrue(self.s.part[0].bonds==bond0 or self.s.part[1].bonds==bond1)
        self.assertTrue(self.s.part[2].bonds==())

        # Check that no additional bonds appear
        self.s.integrator.run(1)
        self.assertTrue(self.s.part[0].bonds==bond0 or self.s.part[1].bonds==bond1)
        self.assertTrue(self.s.part[2].bonds==())


        # Check turning it off
        self.s.collision_detection.set_params(mode=CollisionMode.off)
        self.assertEquals(self.s.collision_detection.mode,CollisionMode.off)
    
    
    @ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_RELATIVE"),"VIRTUAL_SITES not compiled in")
    def test_bind_at_point_of_collision(self):
        self.s.part.clear()
        self.s.part.add(pos=(0,0,0),id=0)
        self.s.part.add(pos=(0.1,0,0),id=1)
        self.s.part.add(pos=(0.1,0.3,0),id=2)

        
        self.s.collision_detection.set_params(mode=CollisionMode.bind_at_point_of_collision,distance=0.11,bond_centers=self.H,bond_vs=self.H2,part_type_vs=1,vs_placement=0.4)
        self.verify_state_after_bind_at_poc()

        # Integrate again and check that nothing has changed
        self.s.integrator.run(1)
        self.verify_state_after_bind_at_poc()

    def verify_state_after_bind_at_poc(self):
        self.s.integrator.run(1,recalc_forces=True)
        bond0=((self.s.bonded_inter[0],1),)
        bond1=((self.s.bonded_inter[0],0),)
        self.assertTrue(self.s.part[0].bonds==bond0 or self.s.part[1].bonds==bond1)
        self.assertTrue(self.s.part[2].bonds==())

        # Check for presence of vs
        vs1=self.s.part[3]
        vs2=self.s.part[4]
        # No additional particles?
        self.assertEquals(self.s.part.highest_particle_id,4)

        # Check for bond betwen vs
        vs_bond1=((self.s.bonded_inter[1],4),)
        vs_bond2=((self.s.bonded_inter[1],3),)
        self.assertTrue(vs1.bonds==vs_bond1 or vs2.bonds==vs_bond2)

        # Vs properties
        self.assertEquals(vs1.virtual,1)
        self.assertEquals(vs2.virtual,1)


        # vs_relative properties
        seen=[]
        for p in vs1,vs2:
          r=p.vs_relative
          rel_to=r[0]
          dist=r[1]
          # Vs is related to one of the particles
          self.assertTrue(rel_to==0 or rel_to==1)
          # The two vs relate to two different particles
          self.assertNotIn(rel_to,seen)
          seen.append(rel_to)

          # Check placement
          if rel_to==0:
            dist_centers=self.s.part[1].pos-self.s.part[0].pos
          else:
            dist_centers=self.s.part[0].pos-self.s.part[1].pos
          expected_pos=self.s.part[rel_to].pos+self.s.collision_detection.vs_placement *dist_centers
          print( p.pos,expected_pos)
          self.assertTrue(np.sqrt(np.sum((p.pos-expected_pos)**2))<=1E-5)
          







        



if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
