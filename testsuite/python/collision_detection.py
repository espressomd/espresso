
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
    
    
    def run_test_bind_at_point_of_collision_for_pos(self,pos):
        self.s.part.clear()
        self.s.part.add(pos=pos+(0,0,0),id=0)
        self.s.part.add(pos=pos+(0.1,0,0),id=1)
        self.s.part.add(pos=pos+(0.1,0.3,0),id=2)

        
        self.s.collision_detection.set_params(mode=CollisionMode.bind_at_point_of_collision,distance=0.11,bond_centers=self.H,bond_vs=self.H2,part_type_vs=1,vs_placement=0.4)
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_state_after_bind_at_poc()


        # Integrate again and check that nothing has changed
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_state_after_bind_at_poc()

    def verify_state_after_bind_at_poc(self):
        bond0=((self.s.bonded_inter[0],1),)
        bond1=((self.s.bonded_inter[0],0),)
        print(bond0,bond1,self.s.part[0].bonds,self.s.part[1].bonds)
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
          print("distance",dist_centers)
          expected_pos=self.s.part[rel_to].pos+self.s.collision_detection.vs_placement *dist_centers
          print("Pos found vs expected:", p.pos,expected_pos)
          self.assertTrue(np.sqrt(np.sum((p.pos-expected_pos)**2))<=1E-5)
    
    @ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_RELATIVE"),"VIRTUAL_SITES not compiled in")
    #@ut.skipIf(s.cell_system.get_state()["n_nodes"]>1,"VS based tests only on a single node")
    def test_bind_at_point_of_collision(self):
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0,0,0)))
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.45,0,0)))

    #@ut.skipIf(not espressomd.has_features("ANGLE_HARMONIC"),"Tests skipped because ANGLE_HARMONIC not compiled in")
    def test_angle_harmonic(self):
        # Setup particles
        self.s.part.clear()
        dx=np.array((1,0,0))
        dy=np.array((0,1,0))
        dz=np.array((0,0,1))
        a=np.array((0.499,0.499,0.499))
        b=a+0.1*dx
        c=a+0.03*dx +0.03*dy
        d=a+0.03*dx -0.03*dy
        e=a-0.1*dx

        self.s.part.add(id=0,pos=a)
        self.s.part.add(id=1,pos=b)
        self.s.part.add(id=2,pos=c)
        self.s.part.add(id=3,pos=d)
        self.s.part.add(id=4,pos=e)


        # Setup bonds
        res=181
        for i in range(0,res,1):
           self.s.bonded_inter[i+2]=Angle_Harmonic(bend=1,phi0=float(i)/(res-1)*np.pi)
        cutoff=0.11
        self.s.collision_detection.set_params(mode=CollisionMode.bind_three_particles,bond_centers=self.H,bond_three_particles=2,three_particle_binding_angle_resolution=res,distance=cutoff)
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_triangle_binding(cutoff,self.s.bonded_inter[2],res)

        # Make sure no extra bonds appear
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_triangle_binding(cutoff,self.s.bonded_inter[2],res)

        # Place the particles in two steps and make sure, the bonds are the same
        self.s.part.clear()
        self.s.part.add(id=0,pos=a)
        self.s.part.add(id=2,pos=c)
        self.s.part.add(id=3,pos=d)
        self.s.integrator.run(0,recalc_forces=True)
        
        self.s.part.add(id=4,pos=e)
        self.s.part.add(id=1,pos=b)
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_triangle_binding(cutoff,self.s.bonded_inter[2],res)

    def verify_triangle_binding(self,distance,first_bond,angle_res):
        # Gather pairs
        n=len(self.s.part)
        angle_res=angle_res-1

        expected_pairs=[]
        for i in range(n):
            for j in range(i+1,n,1):
                if self.s.distance(self.s.part[i],self.s.part[j])<=distance:
                    expected_pairs.append((i,j))
        
        # Find triangles
        # Each elemtn is a particle id, a bond id and two bond partners in ascending order
        expected_angle_bonds=[]
        for i in range(n):
            for j in range(i+1,n,1):
                for k in range(j+1,n,1):
                    # Ref to particles 
                    p_i=self.s.part[i]
                    p_j=self.s.part[j]
                    p_k=self.s.part[k]
                    
                    # Normalized distnace vectors
                    d_ij=p_j.pos-p_i.pos
                    d_ik=p_k.pos-p_i.pos
                    d_jk=p_k.pos-p_j.pos
                    d_ij/=np.sqrt(np.sum(d_ij**2))
                    d_ik/=np.sqrt(np.sum(d_ik**2))
                    d_jk/=np.sqrt(np.sum(d_jk**2))

                    if self.s.distance(p_i,p_j)<=distance and self.s.distance(p_i,p_k)<=distance:
                        id_i=first_bond._bond_id+int(np.round(np.arccos(np.dot(d_ij,d_ik))*angle_res/np.pi))
                        expected_angle_bonds.append((i,id_i,j,k))
                        
                    if self.s.distance(p_i,p_j)<=distance and self.s.distance(p_j,p_k)<=distance:
                        id_j=first_bond._bond_id+int(np.round(np.arccos(np.dot(-d_ij,d_jk))*angle_res/np.pi))
                        expected_angle_bonds.append((j,id_j,i,k))
                    if self.s.distance(p_i,p_k)<=distance and self.s.distance(p_j,p_k)<=distance:
                        id_k=first_bond._bond_id+int(np.round(np.arccos(np.dot(-d_ik,-d_jk))*angle_res/np.pi))
                        expected_angle_bonds.append((k,id_k,i,j))
                       

        # Gather actual pairs and actual triangles
        found_pairs=[]
        found_angle_bonds=[]
        for i in range(n):
            for b in self.s.part[i].bonds:
                if len(b)==2:
                    self.assertEquals(b[0]._bond_id,self.H._bond_id)
                    found_pairs.append(tuple(sorted((i,b[1]))))
                elif len(b)==3:
                    partners=sorted(b[1:])
                    found_angle_bonds.append((i,b[0]._bond_id,partners[0],partners[1]))
                else:
                    raise Exception("There should be only 2 and three particle bonds")
        
        # The order between expected and found bonds does not malways match
        # because collisions occur in random order. Sort stuff
        found_pairs=sorted(found_pairs)
        found_angle_bonds=sorted(found_angle_bonds)
        expected_angle_bonds=sorted(expected_angle_bonds)
        self.assertEquals(expected_pairs,found_pairs)
        
        if not  expected_angle_bonds == found_angle_bonds:
            # Verbose info
            print("expected:",expected_angle_bonds)
            missing=[]
            for b in expected_angle_bonds:
                if b in found_angle_bonds:
                    found_angle_bonds.remove(b)
                else:
                    missing.append(b)
            print("missing",missing)
            print("extra:",found_angle_bonds)
            print()
        
        self.assertEquals(expected_angle_bonds,found_angle_bonds)
            
                
                    

                  


            
        
        








        



if __name__ == "__main__":
    ut.main()
