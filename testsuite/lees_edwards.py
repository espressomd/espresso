#!/usr/bin/env python

from __future__ import print_function

import espressomd as md
import unittest as ut

@ut.skipIf(not md.has_features(['LEES_EDWARDS']),
'Features not available, skipping test!')

class Lees_Edwards_test(ut.TestCase):

  def test(self):
    
    #Systemclass
    system = md.System()
    
    #allowed deviation from analytical results
    tol = 10e-15
    
    #Simulation box and integration parameters
    box_l = 5.0
    system.box_l = [box_l, box_l, box_l]
    time_step = 1.0
    system.time_step = time_step 
    system.cell_system.skin = 0.4
    
    #Check if coordinates are folded correctly
    #Particle placement 
    system.lees_edwards_offset = 0.1
    system.part.add(pos=[2.5, 6.0, 2.5], v=[0.05, 0.1, 0.05], id=0, type=0)
    pos_difference_x = system.part[0].pos_folded[0] - 2.4
    pos_difference_y = system.part[0].pos_folded[1] - 1.0
    pos_difference_z = system.part[0].pos_folded[2] - 2.5
    vel_difference_x = system.part[0].v[0] - 0.05
    vel_difference_y = system.part[0].v[1] - 0.1
    vel_difference_z = system.part[0].v[2] - 0.05
    
    self.assertTrue(pos_difference_x < tol)
    self.assertTrue(pos_difference_y < tol)
    self.assertTrue(pos_difference_z < tol)
    self.assertTrue(vel_difference_x < tol)
    self.assertTrue(vel_difference_y < tol)
    self.assertTrue(vel_difference_z < tol)

    print(pos_difference_x, pos_difference_y, pos_difference_z, vel_difference_x, vel_difference_y, vel_difference_z)
    
    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    #Check position offset in x-direction while crossing upper boundary
    #Particle placement
    system.part.add(pos=[2.5, 4.8, 2.5], v=[0.0, 0.1, 0.0], id=0, type=0)
    
    #Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[1] - 0.1
      print(pos_difference, vel_difference)

      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)

    for i in range(3):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.4
      vel_difference = system.part[0].v[1] -0.1
      print(pos_difference, vel_difference)
   
      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)

    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    #Check position offset in x-direction while crossing lower boundary
    #Particle placement
    system.part.add(pos=[2.5, 0.2, 2.5], v=[0.0, -0.1, 0.0], id=0, type=0)
    
    #Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[1] - -0.1
      print(pos_difference, vel_difference)
      
      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)
    
    for i in range(3):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.6
      vel_difference = system.part[0].v[1] - -0.1
      print(pos_difference, vel_difference)

      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)

    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    #Check velocity offset with constant shear rate while crossing upper boundary
    #Particle placement
    system.part.add(pos=[2.5, 4.8, 2.5], v=[0.0, 0.1, 0.0], id=0, type=0)
    
    #Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[0] + 0.0
      #print(system.part[0].pos_folded, system.part[0].v)
      print(pos_difference, vel_difference)
      
      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)
    
    for i in range(3,6):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.6 + i*0.1
      vel_difference = system.part[0].v[0] + 0.1
      #print(system.part[0].pos_folded, system.part[0].v)
      print(pos_difference, vel_difference)
      
      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)
    
    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    #Check position offset in x-direction while crossing lower boundary
    #Particle placement
    system.part.add(pos=[2.5, 0.2, 2.5], v=[0.0, -0.1, 0.0], id=0, type=0)
    
    #Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[0] - 0.0
      #print(system.part[0].pos_folded, system.part[0].v)
      print(pos_difference, vel_difference)
      
      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)
    
    for i in range(3,6):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.4 - i*0.1
      vel_difference = system.part[0].v[0] - 0.1
      #print(system.part[0].pos_folded, system.part[0].v)
      print(pos_difference, vel_difference)
      
      self.assertTrue(pos_difference < tol)
      self.assertTrue(vel_difference < tol)
    
    system.part.clear()
    system.lees_edwards_offset = 0.0
    
if __name__ == "__main__":
  ut.main()
