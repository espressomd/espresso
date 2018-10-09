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

        # Protocol should be off by default
        self.assertEqual(system.lees_edwards[0], "off")

        # One shouldn't be able to set a random protocol
        with self.assertRaises(Exception):
            system.lees_edwards = ["oscillatory_shear", 50.]

        # Check if the setted protocol is stored correctly
        system.lees_edwards = ["steady_shear", 1.2, 2, 0]
        self.assertEqual(system.lees_edwards[0], "steady_shear")
        self.assertEqual(system.lees_edwards[1], 1.2)
        self.assertEqual(system.lees_edwards[3], 2)
        self.assertEqual(system.lees_edwards[4], 0)

        system.lees_edwards = ["oscillatory_shear", 1.2, 5.6, 0, 1]
        self.assertEqual(system.lees_edwards[0], "oscillatory_shear")
        self.assertEqual(system.lees_edwards[1], 1.2)
        self.assertEqual(system.lees_edwards[2], 5.6)
        self.assertEqual(system.lees_edwards[5], 0)
        self.assertEqual(system.lees_edwards[6], 1)

        system.time = 0.0

        # Check if the offset is determined correctly
        frequency = 3.7
        omega = 2 * np.pi * frequency
        amplitude = 1.6
        sheardir = 0 
        shearplanenormal = 1
        system.lees_edwards = ["oscillatory_shear", omega, amplitude, sheardir, shearplanenormal]
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

        velocity = 0.1
        dir = [0, 1, 2]

        for sheardir in dir:
          for shearplanenormal in dir:
            if sheardir != shearplanenormal:

              system.time = 0.0
              system.lees_edwards = ["steady_shear", velocity, sheardir, shearplanenormal]

              pos = np.full([3], 2.5)
              pos[shearplanenormal] = system.box_l[shearplanenormal] - 0.1
              vel = np.zeros([3])
              vel[shearplanenormal] = 0.1
              system.part.add(pos=pos, v=vel, id=0, type=0)

              pos_change = np.zeros([3])
              pos_change[sheardir] = -0.5*system.time_step*velocity
              pos_change[shearplanenormal] = velocity*system.time_step
              vel_change = np.zeros([3])
              vel_change[sheardir] = -velocity

              expected_pos = system.part[0].pos + pos_change
              expected_vel = system.part[0].v + vel_change

              system.integrator.run(1)

              np.testing.assert_almost_equal(system.part[0].pos, expected_pos)
              np.testing.assert_almost_equal(system.part[0].v, expected_vel)

              system.part.clear()

              # Test for lower boundary

              system.time = 0.0
              pos = np.full([3], 2.5)
              pos[shearplanenormal] = 0.1
              vel = np.zeros([3])
              vel[shearplanenormal] = -0.1
              system.part.add(pos=pos, v=vel, id=0, type=0)

              pos_change = np.zeros([3])
              pos_change[sheardir] = 0.5*system.time_step*velocity
              pos_change[shearplanenormal] = -velocity*system.time_step
              vel_change = np.zeros([3])
              vel_change[sheardir] = velocity

              expected_pos = system.part[0].pos + pos_change
              expected_vel = system.part[0].v + vel_change

              system.integrator.run(1)

              np.testing.assert_almost_equal(system.part[0].pos, expected_pos)
              np.testing.assert_almost_equal(system.part[0].v, expected_vel)

              system.part.clear()

if __name__ == "__main__":
    ut.main()
