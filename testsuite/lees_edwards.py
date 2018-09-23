from __future__ import print_function

import espressomd
import unittest as ut
import numpy as np

@ut.skipIf(not espressomd.has_features(["LEES_EDWARDS"]),
  "Feature not available, skipping test!")
class LeesEdwards(ut.TestCase):

  system = espressomd.System(box_l=[5.0, 5.0, 5.0])
  system.cell_system.skin = 0.0
  system.cell_system.set_n_square()
  system.set_random_state_PRNG()

  tol = 1e-15
  time_step = 1.0
  system.time_step = time_step

  def test_Protocol(self):
    """This test calculates if the offset and velocity of the Lees Edwards
    function are calculated correctly based on the input variables"""

    system = self.system

    #Protocol should be off by default
    self.assertEqual(system.lees_edwards[0], "off")

    #One shouldn't be able to set a random protocol
    with self.assertRaises(Exception):
      system.lees_edwards = ["oscillatory_shear", 50.]

    #Check if the setted protocol is stored correctly
    system.lees_edwards = ["steady_shear", 1.2]
    self.assertEqual(system.lees_edwards[0], "steady_shear")
    self.assertEqual(system.lees_edwards[1], 1.2)

    system.lees_edwards = ["oscillatory_shear", 1.2, 5.6]
    self.assertEqual(system.lees_edwards[0], "oscillatory_shear")
    self.assertEqual(system.lees_edwards[1], 1.2)
    self.assertEqual(system.lees_edwards[2], 5.6)

    system.time = 5.0

    #Check if the offset is determined correctly
    frequency = 3.7
    omega = 2*np.pi*frequency
    amplitude = 1.6
    system.lees_edwards = ["oscillatory_shear", frequency, amplitude]
    self.assertEqual(system.lees_edwards[0], "oscillatory_shear")
    for time in np.arange(0.0, 100.0, 10.0):
      offset = amplitude * np.sin(frequency * time)
      velocity = omega * amplitude * np.cos(omega * time)

      np.testing.assert_equal(system.lees_edwards[3], velocity)
      np.testing.assert_equal(system.lees_edwards[4], offset)
      system.integrator.run(10)

  def BoundaryCrossingTest(self):
    """The test calculates a particle's position for different lees_edwards_offsets
    while crossing the upper and lower boundary in the y-direction. The
    obtained positions as well as the velocity of the particle is checked via a 
    comparison with the expected value."""

    # Check if coordinates are folded correctly
    system.lees_edwards = ["steady_shear", 0.1]
    
    # Particle placement
    system.part.add(pos=[2.5, 6.0, 2.5], v=[0.05, 0.1, 0.05], id=0, type=0)
    
    pos_difference_x = system.part[0].pos_folded[0] - 2.4
    pos_difference_y = system.part[0].pos_folded[1] - 1.0
    pos_difference_z = system.part[0].pos_folded[2] - 2.5
    vel_difference_x = system.part[0].v[0] - 0.05
    vel_difference_y = system.part[0].v[1] - 0.1
    vel_difference_z = system.part[0].v[2] - 0.05
    
    self.assertLess(pos_difference_x, tol)
    self.assertLess(pos_difference_y, tol)
    self.assertLess(pos_difference_z, tol)
    self.assertLess(vel_difference_x, tol)
    self.assertLess(vel_difference_y, tol)
    self.assertLess(vel_difference_z, tol)

    system.part.clear()
    system.lees_edwards_offset = 0.0

    # Check position offset in x-direction while crossing upper boundary
    # Particle placement
    system.part.add(pos=[2.5, 4.8, 2.5], v=[0.0, 0.1, 0.0], id=0, type=0)
    
    # Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[1] - 0.1

      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)

    for i in range(3):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.4
      vel_difference = system.part[0].v[1] - 0.1
   
      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)

    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    # Check position offset in x-direction while crossing lower boundary
    # Particle placement
    system.part.add(pos=[2.5, 0.2, 2.5], v=[0.0, -0.1, 0.0], id=0, type=0)
    
    # Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[1] - -0.1
      
      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)
    
    for i in range(3):
      system.integrator.run(1)
      system.lees_edwards_offset = 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.6
      vel_difference = system.part[0].v[1] - -0.1

      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)

    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    # Check velocity offset with constant shear rate while crossing upper
    # boundary
    # Particle placement
    system.part.add(pos=[2.5, 4.8, 2.5], v=[0.0, 0.1, 0.0], id=0, type=0)
    
    # Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[0] + 0.0
      
      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)
    
    for i in range(3, 6):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.6 + i * 0.1
      vel_difference = system.part[0].v[0] + 0.1
      
      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)
    
    system.part.clear()
    system.lees_edwards_offset = 0.0
    
    # Check position offset in x-direction while crossing lower boundary
    # Particle placement
    system.part.add(pos=[2.5, 0.2, 2.5], v=[0.0, -0.1, 0.0], id=0, type=0)
    
    # Integration
    for i in range(2):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.5
      vel_difference = system.part[0].v[0] - 0.0
      
      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)
    
    for i in range(3, 6):
      system.integrator.run(1)
      system.lees_edwards_offset += 0.1
      pos_difference = system.part[0].pos_folded[0] - 2.4 - i * 0.1
      vel_difference = system.part[0].v[0] - 0.1
      
      self.assertLess(pos_difference, tol)
      self.assertLess(vel_difference, tol)
    
    system.part.clear()
    system.lees_edwards_offset = 0.0
    
if __name__ == "__main__":
  ut.main()
