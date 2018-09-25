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

  def test_a_Protocol(self):
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

    system.time = 0.0

    #Check if the offset is determined correctly
    frequency = 3.7
    omega = 2*np.pi*frequency
    amplitude = 1.6
    system.lees_edwards = ["oscillatory_shear", omega, amplitude]
    self.assertEqual(system.lees_edwards[0], "oscillatory_shear")
    for time in np.arange(10., 100.0, 10.0):
      system.integrator.run(10)
      offset = amplitude * np.sin(omega * time)
      velocity = omega * amplitude * np.cos(omega * time)

      np.testing.assert_almost_equal(system.lees_edwards[3], velocity)
      np.testing.assert_almost_equal(system.lees_edwards[4], offset)

  def test_b_BoundaryCrossing(self):
    """A particle crosses the upper and lower boundary to test if position
    and velocity are updated correctly."""

    system = self.system
    tol = self.tol
    system.time = 0.0

    # Set up a one particle system and check the position offset after crossing the boundary
    # Test for upper boundary
    system.lees_edwards = ["steady_shear", 0.1]
    system.part.add(pos=[2.5, 4.9, 2.5], v=[0.0, 0.1, 0.0], id=0, type=0)

    system.integrator.run(1)

    expected_pos = [2.45, 5.0, 2.5]
    expected_vel = [-0.1, 0.1, 0.0]
    
    np.testing.assert_almost_equal(system.part[0].pos, expected_pos)
    np.testing.assert_almost_equal(system.part[0].v, expected_vel)

    system.part.clear()

    # Test for lower boundary
    
    system.time = 0.0
    system.part.add(pos=[2.5, 0.1, 2.5], v=[0.0, -0.1, 0.0], id=0, type=0)

    system.integrator.run(1)

    expected_pos = [2.55, 0., 2.5]
    expected_vel = [0.1, -0.1, 0.0]
    
    print(system.part[0].pos_folded)

    np.testing.assert_almost_equal(system.part[0].pos, expected_pos)
    np.testing.assert_almost_equal(system.part[0].v, expected_vel)

if __name__ == "__main__":
  ut.main()
