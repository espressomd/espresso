
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
from espressomd.interactions import HarmonicBond,AngleHarmonic
import numpy as np
from random import shuffle

@ut.skipIf(not espressomd.has_features("COLLISION_DETECTION"),"Required features not compiled in")
class CollisionDetection(ut.TestCase):
    """Tests interface and functionality of the collision detection / dynamic binding"""

    s = espressomd.System(box_l = [10.0, 10.0, 10.0])
    s.seed  = s.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(seed=s.seed)
    if espressomd.has_features("VIRTUAL_SITES"): 
        from espressomd.virtual_sites import VirtualSitesRelative
        s.virtual_sites=VirtualSitesRelative()

    H = HarmonicBond(k=5000,r_0=0.1)
    H2 = HarmonicBond(k=25000,r_0=0.02)
    s.bonded_inter.add(H)
    s.bonded_inter.add(H2)
    s.time_step=0.001
    s.cell_system.skin=0.05
    s.min_global_cut=0.2

    part_type_to_attach_vs_to=0
    part_type_vs=1
    part_type_to_be_glued=2
    part_type_after_glueing=3 
    other_type=5

    def test_00_interface_and_defaults(self):
        # Is it off by default
        self.assertEqual(self.s.collision_detection.mode,"off")
        # Make sure params cannot be set individually
        with self.assertRaises(Exception):
            self.s.collision_detection.mode="bind_centers" 
        
        # Verify exception throwing for unknown collision modes
        with self.assertRaises(Exception):
            self.s.collision_detection.set_params(mode=0)
            self.s.collision_detection.set_params(mode="blahblah")
        
        # That should work
        self.s.collision_detection.set_params(mode="off")
        self.assertEqual(self.s.collision_detection.mode,"off")

    def test_bind_centers(self):
        # Check that it leaves particles alone, wehn off
        self.s.collision_detection.set_params(mode="off")
        
        self.s.part.clear()
        self.s.part.add(pos=(0,0,0),id=0)
        self.s.part.add(pos=(0.1,0,0),id=1)
        self.s.part.add(pos=(0.1,0.3,0),id=2)
        self.s.integrator.run(0)
        self.assertEqual(self.s.part[0].bonds,())
        self.assertEqual(self.s.part[1].bonds,())
        self.assertEqual(self.s.part[2].bonds,())

        # Check that it cannot be activated 
        self.s.collision_detection.set_params(mode="bind_centers",distance=0.11,bond_centers=self.H)
        self.s.integrator.run(1,recalc_forces=True)
        bond0=((self.s.bonded_inter[0],1),)
        bond1=((self.s.bonded_inter[0],0),)
        self.assertTrue(self.s.part[0].bonds==bond0 or self.s.part[1].bonds==bond1)
        self.assertEqual(self.s.part[2].bonds,())

        # Check that no additional bonds appear
        self.s.integrator.run(1)
        self.assertTrue(self.s.part[0].bonds==bond0 or self.s.part[1].bonds==bond1)
        self.assertEqual(self.s.part[2].bonds,())


        # Check turning it off
        self.s.collision_detection.set_params(mode="off")
        self.assertEqual(self.s.collision_detection.mode,"off")
    
    
    def run_test_bind_at_point_of_collision_for_pos(self,*positions):
        positions=list(positions)
        shuffle(positions)
        self.s.part.clear()
        # Place particle which should not take part in collisions
        p=self.s.part.add(pos=(0.1,0.3,0))
        for pos in positions:
            p1=self.s.part.add(pos=pos+(0,0,0))
            p2=self.s.part.add(pos=pos+(0.1,0,0))
            if self.s.distance(p1,p) <0.12 or self.s.distance(p2,p)<0.12:
                raise Exception("Test particle too close to particle, which should not take part in collision")


        
        # 2 non-virtual + 2 virtual + one that doesn't tkae part
        expected_np=4*len(positions)+1

        self.s.collision_detection.set_params(mode="bind_at_point_of_collision",distance=0.11,bond_centers=self.H,bond_vs=self.H2,part_type_vs=1,vs_placement=0.4)
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)


        # Integrate again and check that nothing has changed
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Check that nothing explodes, when the particles are moved.
        # In particular for parallel simulations
        self.s.thermostat.set_langevin(kT=0,gamma=0.01)
        self.s.part[:].v=0.05,0.01,0.15
        self.s.integrator.run(3000)
        self.verify_state_after_bind_at_poc(expected_np)


    def verify_state_after_bind_at_poc(self,expected_np):
        self.assertEqual(len(self.s.part),expected_np)

        # At the end of test, this list should be empty
        parts_not_accounted_for=list(range(expected_np))
        
        # Collect pairs of non-virtual-particles found
        non_virtual_pairs=[]

        # We traverse particles. We look for a vs with a bond to find the other vs.
        # From the two vs we find the two non-virtual particles
        for p in self.s.part:
            # Skip non-virtual
            if p.virtual==0:
                continue
            # Skip vs that doesn't have a bond
            if p.bonds==():
                continue
            # Parse the bond
            self.assertEqual(len(p.bonds),1)
            # Bond type
            self.assertEqual(p.bonds[0][0],self.H2)
            # get partner
            p2 =self.s.part[p.bonds[0][1]]
            # Is that really a vs
            self.assertEqual(p2.virtual,1)
            # Get base particles
            base_p1=self.s.part[p.vs_relative[0]]
            base_p2=self.s.part[p2.vs_relative[0]]
            # Take note of accounted-for particles
            for _p in p,p2,base_p1,base_p2:
                parts_not_accounted_for.remove(_p.id)
            self.verify_bind_at_poc_pair(base_p1,base_p2,p,p2)
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for),1)
        p=self.s.part[parts_not_accounted_for[0]]
        self.assertEqual(p.virtual,0)
        self.assertEqual(p.bonds,())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for,[])


    def verify_bind_at_poc_pair(self,p1,p2,vs1,vs2):
        bond_p1=((self.s.bonded_inter[0],p2.id),)
        bond_p2=((self.s.bonded_inter[0],p1.id),)
        self.assertTrue(p1.bonds==bond_p1 or p2.bonds==bond_p2)

        # Check for presence of vs
        # Check for bond betwen vs
        bond_vs1=((self.s.bonded_inter[1],vs2.id),)
        bond_vs2=((self.s.bonded_inter[1],vs1.id),)
        self.assertTrue(vs1.bonds==bond_vs1 or vs2.bonds==bond_vs2)

        # Vs properties
        self.assertEqual(vs1.virtual,1)
        self.assertEqual(vs2.virtual,1)


        # vs_relative properties
        seen=[]
        for p in vs1,vs2:
          r=p.vs_relative
          rel_to=r[0]
          dist=r[1]
          # Vs is related to one of the particles
          self.assertTrue(rel_to==p1.id or rel_to==p2.id)
          # The two vs relate to two different particles
          self.assertNotIn(rel_to,seen)
          seen.append(rel_to)

          # Check placement
          if rel_to==p1.id:
            dist_centers=np.copy(p2.pos-p1.pos)
          else:
            dist_centers=p1.pos-p2.pos
          expected_pos=self.s.part[rel_to].pos_folded+self.s.collision_detection.vs_placement *dist_centers
          np.testing.assert_allclose(np.copy(p.pos_folded),expected_pos,atol=1E-4)
    
    @ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_RELATIVE"),"VIRTUAL_SITES not compiled in")
    #@ut.skipIf(s.cell_system.get_state()["n_nodes"]>1,"VS based tests only on a single node")
    def test_bind_at_point_of_collision(self):
        # Single collision head node
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0,0,0)))
        # Single collision, mixed
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.45,0,0)))
        # Single collision, non-head-node
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.7,0,0)))

        # head-node + mixed
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0,0,0)),np.array((0.45,0,0)))
        # Mixed + other node
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.45,0,0)),np.array((0.7,0,0)))
        # Head + other
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.0,0,0)),np.array((0.7,0,0)))
        # Head + mixed + other
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.2,0,0)),np.array((0.95,0,0)),np.array((0.7,0,0)))
    
    
    
    def run_test_glue_to_surface_for_pos(self,*positions):
        positions=list(positions)
        shuffle(positions)
        self.s.part.clear()
        # Place particle which should not take part in collisions
        # In this case, it is skipped, because it is of the wrong type,
        # even if it is within range for a collision
        p=self.s.part.add(pos=positions[0],type=self.other_type)
        for pos in positions:
            # Since this is non-symmetric, we randomize order
            if np.random.random()>.5:
                p1=self.s.part.add(pos=pos+(0,0,0),type=self.part_type_to_attach_vs_to)
                p2=self.s.part.add(pos=pos+(0.1,0,0),type=self.part_type_to_be_glued)
            else:
                p2=self.s.part.add(pos=pos+(0.1,0,0),type=self.part_type_to_be_glued)
                p1=self.s.part.add(pos=pos+(0,0,0),type=self.part_type_to_attach_vs_to)

        
        # 2 non-virtual + 1 virtual + one that doesn't takekae part
        expected_np=3*len(positions)+1

        self.s.collision_detection.set_params(
          mode="glue_to_surface",distance=0.11,distance_glued_particle_to_vs=0.02, bond_centers=self.H,bond_vs=self.H2,part_type_vs=self.part_type_vs,part_type_to_attach_vs_to=self.part_type_to_attach_vs_to,part_type_to_be_glued=self.part_type_to_be_glued,part_type_after_glueing=self.part_type_after_glueing)
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_state_after_glue_to_surface(expected_np)


        # Integrate again and check that nothing has changed
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_state_after_glue_to_surface(expected_np)

        # Check that nothing explodes, when the particles are moved.
        # In particular for parallel simulations
        self.s.thermostat.set_langevin(kT=0,gamma=0.01)
        self.s.part[:].v=0.05,0.01,0.15
        self.s.integrator.run(3000)
        self.verify_state_after_glue_to_surface(expected_np)


    def verify_state_after_glue_to_surface(self,expected_np):
        self.assertEqual(len(self.s.part),expected_np)

        # At the end of test, this list should be empty
        parts_not_accounted_for=list(range(expected_np))
        
        # We traverse particles. We look for a vs, get base particle from there
        # and prtner particle via bonds
        for p in self.s.part:
            # Skip non-virtual
            if p.virtual==0:
                continue
            # The vs shouldn't have bonds
            self.assertEqual(p.bonds,())

            # Get base particles
            base_p=self.s.part[p.vs_relative[0]]
            
            # Get bound particle
            # There is a bond between the base particle and the bound particle
            # but we have no guarantee, on where its stored
            # 1. On the base particle of the vs
            p2=None
            if len(base_p.bonds)==1:
                self.assertEqual(base_p.bonds[0][0],self.H)
                p2=self.s.part[base_p.bonds[0][1]]
            else:
                # We need to go through all particles to find it
                for candidate in self.s.part:
                    if candidate.id not in parts_not_accounted_for:
                        continue
                    if len(candidate.bonds)>=1:
                        for b in candidate.bonds:
                            if b[0] == self.H and b[1]==base_p.id:
                                p2=candidate
                if p2==None:
                    raise Exception("Bound particle not found")
            # Take note of accounted-for particles
            parts_not_accounted_for.remove(base_p.id)
            parts_not_accounted_for.remove(p.id)
            parts_not_accounted_for.remove(p2.id)
            self.verify_glue_to_surface_pair(base_p,p,p2)
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for),1)
        p=self.s.part[parts_not_accounted_for[0]]
        self.assertEqual(p.virtual,0)
        self.assertEqual(p.type,self.other_type)
        self.assertEqual(p.bonds,())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for,[])


    def verify_glue_to_surface_pair(self,base_p,vs,bound_p):
        # Check all types
        self.assertEqual(base_p.type,self.part_type_to_attach_vs_to)
        self.assertEqual(vs.type,self.part_type_vs)
        self.assertEqual(bound_p.type,self.part_type_after_glueing)
        
        # Bound particle should have a bond to vs. It can additionally have a bond
        # to the base particle
        bond_to_vs_found=False
        for b in bound_p.bonds: 
            if b[0] == self.H2: 
                # bond to vs
                self.assertEqual(b,(self.H2,vs.id))
                bond_to_vs_found=True
        self.assertEqual(bond_to_vs_found,True)
        # Vs should not have a bond
        self.assertEqual(vs.bonds,())
        
        # Vs properties
        self.assertEqual(vs.virtual,1)
        self.assertEqual(vs.vs_relative[0],base_p.id)
        
        # Distance vs,bound_p
        self.assertAlmostEqual(self.s.distance(vs,bound_p),0.02,places=3)
        self.assertAlmostEqual(self.s.distance(base_p,bound_p),0.1,places=3)
        self.assertAlmostEqual(self.s.distance(base_p,vs),0.08,places=3)

        # base_p,vs,bound_p on a line
        self.assertGreater(np.dot(self.s.distance_vec(base_p,vs),self.s.distance_vec(base_p,bound_p))/self.s.distance(base_p,vs)/self.s.distance(base_p,bound_p),0.99)


       


    @ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_RELATIVE"),"Skipped due to missing VIRTUAL_SITES_RELATIVE")
    def test_glue_to_surface(self):
        # Single collision head node
        self.run_test_glue_to_surface_for_pos(np.array((0,0,0)))
        # Single collision, mixed
        self.run_test_glue_to_surface_for_pos(np.array((0.45,0,0)))
        # Single collision, non-head-node
        self.run_test_glue_to_surface_for_pos(np.array((0.7,0,0)))

        # head-node + mixed
        self.run_test_glue_to_surface_for_pos(np.array((0,0,0)),np.array((0.45,0,0)))
        # Mixed + other node
        self.run_test_glue_to_surface_for_pos(np.array((0.45,0,0)),np.array((0.7,0,0)))
        # Head + other
        self.run_test_glue_to_surface_for_pos(np.array((0.0,0,0)),np.array((0.7,0,0)))
        # Head + mixed + other
        self.run_test_glue_to_surface_for_pos(np.array((0.2,0,0)),np.array((0.95,0,0)),np.array((0.7,0,0)))





    

    #@ut.skipIf(not espressomd.has_features("AngleHarmonic"),"Tests skipped because AngleHarmonic not compiled in")
    def test_AngleHarmonic(self):
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
           self.s.bonded_inter[i+2]=AngleHarmonic(bend=1,phi0=float(i)/(res-1)*np.pi)
        cutoff=0.11
        self.s.collision_detection.set_params(mode="bind_three_particles",bond_centers=self.H,bond_three_particles=2,three_particle_binding_angle_resolution=res,distance=cutoff)
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
        self.s.cell_system.set_domain_decomposition()
        self.s.integrator.run(0,recalc_forces=True)
        self.verify_triangle_binding(cutoff,self.s.bonded_inter[2],res)
        self.s.cell_system.set_n_square()
        self.s.part[:].bonds=()
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
                    d_ij=np.copy(p_j.pos-p_i.pos)
                    d_ik=np.copy(p_k.pos-p_i.pos)
                    d_jk=np.copy(p_k.pos-p_j.pos)
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
                    self.assertEqual(b[0]._bond_id,self.H._bond_id)
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
        self.assertEqual(expected_pairs,found_pairs)
        
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
        
        self.assertEqual(expected_angle_bonds,found_angle_bonds)
            
                  









if __name__ == "__main__":
    ut.main()
