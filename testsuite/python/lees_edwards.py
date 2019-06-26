from __future__ import print_function

import espressomd
from espressomd.interactions import HarmonicBond

import unittest as ut
import numpy as np


@ut.skipIf(not espressomd.has_features(["LEES_EDWARDS"]),
           "Feature not available, skipping test!")
class LeesEdwards(ut.TestCase):

    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists = False)
    system.set_random_state_PRNG()

    time_step = 0.5
    system.time_step = time_step

    def test_a_Protocol(self):
        """This test calculates if the offset and velocity of the Lees Edwards
        function are calculated correctly based on the input variables"""

        system = self.system

        # Protocol should be off by default
        self.assertEqual(self.system.lees_edwards.type, "off")

        # One shouldn't be able to set a random protocol
        with self.assertRaises(Exception):
            system.lees_edwards.set_params(
                type="oscillatory_shear", velocity=50)

        # Check if the setted protocol is stored correctly
        system.lees_edwards.set_params(
            type="steady_shear", velocity=1.2, sheardir=2, shearplanenormal=0)
        self.assertEqual(system.lees_edwards.type, "steady_shear")
        self.assertAlmostEqual(system.lees_edwards.velocity, 1.2)
        self.assertEqual(system.lees_edwards.sheardir, 2)
        self.assertEqual(system.lees_edwards.shearplanenormal, 0)

        system.lees_edwards.set_params(
            type="oscillatory_shear", frequency=1.2, amplitude=5.5, sheardir=0, shearplanenormal=1)
        self.assertEqual(system.lees_edwards.type, "oscillatory_shear")
        self.assertAlmostEqual(system.lees_edwards.frequency, 1.2)
        self.assertAlmostEqual(system.lees_edwards.amplitude, 5.5)
        self.assertEqual(system.lees_edwards.sheardir, 0)
        self.assertEqual(system.lees_edwards.shearplanenormal, 1)

        system.time = 0.0

        # Check if the offset is determined correctly
        frequency = np.sqrt(2.0)
        omega = 2 * np.pi * frequency
        amplitude = 1.3
        sheardir = 0
        shearplanenormal = 1
        system.lees_edwards.set_params(type="oscillatory_shear", frequency=omega,
                                       amplitude=amplitude, sheardir=sheardir, shearplanenormal=shearplanenormal)
        self.assertEqual(system.lees_edwards.type, "oscillatory_shear")
        for time in np.arange(system.time_step*10.0, 100.0, system.time_step*10):
            system.integrator.run(10)
            offset = amplitude * np.sin(omega * time)
            velocity = omega * amplitude * np.cos(omega * time)

            print(system.lees_edwards.velocity, velocity)

            self.assertAlmostEqual(system.lees_edwards.velocity, velocity)
            self.assertAlmostEqual(system.lees_edwards.offset, offset)

    def test_b_BoundaryCrossing(self):
        """A particle crosses the upper and lower boundary to test if position
        and velocity are updated correctly."""

        system = self.system
        system.part.clear()
        system.time = 0.0

        # Set up a one particle system and check the position offset after crossing the boundary
        # Test for upper boundary

        velocity = 0.1
        dir = [0, 1, 2]

        for sheardir in dir:
            for shearplanenormal in dir:
                if sheardir != shearplanenormal:

                    system.time = 0.0
                    system.lees_edwards.set_params(
                        type="steady_shear", velocity=velocity, sheardir=sheardir, shearplanenormal=shearplanenormal)

                    pos = np.full([3], 2.5)
                    pos[shearplanenormal] = system.box_l[shearplanenormal] - 0.05
                    vel = np.zeros([3])
                    vel[shearplanenormal] = 0.2
                    system.part.add(pos=pos, v=vel, id=0, type=0)

                    pos_change = np.zeros([3])
                    pos_change[sheardir] = -system.time_step*0.5*velocity
                    pos_change[shearplanenormal] = vel[shearplanenormal]*system.time_step
                    vel_change = np.zeros([3])
                    vel_change[sheardir] = -velocity

                    expected_pos = system.part[0].pos + pos_change
                    expected_vel = system.part[0].v + vel_change

                    system.integrator.run(1)

                    np.testing.assert_almost_equal(
                        system.part[0].pos, expected_pos)
                    np.testing.assert_almost_equal(
                        system.part[0].v, expected_vel)

                    system.part.clear()

                    # Test for lower boundary

                    system.time = 0.0
                    pos = np.full([3], 2.5)
                    pos[shearplanenormal] = 0.05
                    vel = np.zeros([3])
                    vel[shearplanenormal] = -0.2
                    system.part.add(pos=pos, v=vel, id=0, type=0)

                    pos_change = np.zeros([3])
                    pos_change[sheardir] = system.time_step*0.5*velocity
                    pos_change[shearplanenormal] = vel[shearplanenormal]*system.time_step
                    vel_change = np.zeros([3])
                    vel_change[sheardir] = velocity

                    expected_pos = system.part[0].pos + pos_change
                    expected_vel = system.part[0].v + vel_change

                    system.integrator.run(1)

                    np.testing.assert_almost_equal(
                        system.part[0].pos, expected_pos)
                    np.testing.assert_almost_equal(
                        system.part[0].v, expected_vel)

                    system.part.clear()

    def test_c_interactions(self):
        """We place two particles crossing a boundary and connect them with an unbonded
           and bonded interaction. We test with the resulting stress tensor if the offset
           is included properly"""

        system = self.system
        system.part.clear()

        dir = [0, 1, 2]

        for sheardir in dir:
            for shearplanenormal in dir:
                if sheardir != shearplanenormal:

                    offset = 1.0
                    system.lees_edwards.set_params(
                        type='step', offset=offset, sheardir=sheardir, shearplanenormal=shearplanenormal)

                    pos1 = np.full([3], 2.5)
                    pos1[shearplanenormal] = 4.75
                    pos2 = np.full([3], 2.5)
                    pos2[shearplanenormal] = 0.25

                    system.part.add(id=0, pos=pos1, fix=[1, 1, 1])
                    system.part.add(id=1, pos=pos2, fix=[1, 1, 1])

                    r = system.part[1].pos - system.part[0].pos
                    r[sheardir] += offset
                    r[shearplanenormal] += system.box_l[0]

                    k = 1000
                    r_cut = 1.5
                    harm = HarmonicBond(k=k, r_0=0.0)
                    system.bonded_inter.add(harm)
                    system.part[0].add_bond((harm, 1))
                    system.non_bonded_inter[0, 0].soft_sphere.set_params(
                        a=k, n=-2, cutoff=r_cut)

                    system.integrator.run(0, recalc_forces=True)

                    simulated_bondedstress = np.abs(
                        system.analysis.stress_tensor()['bonded'])
                    simulated_nonbondedstress = np.abs(
                        system.analysis.stress_tensor()['non_bonded'])

                    analytical_bondedstress = np.zeros([3, 3])
                    for i in range(3):
                        for j in range(3):
                            analytical_bondedstress[i, j] += k*r[i]*r[j]
                    analytical_bondedstress /= (system.box_l[0]**3.0)

                    analytical_nonbondedstress = np.zeros([3, 3])
                    for i in range(3):
                        for j in range(3):
                            analytical_nonbondedstress[i, j] += 2*k*r[i]*r[j]
                    analytical_nonbondedstress /= (system.box_l[0]**3.0)

                    np.testing.assert_array_equal(
                        simulated_bondedstress, analytical_bondedstress)
                    np.testing.assert_array_equal(
                        simulated_nonbondedstress, analytical_nonbondedstress)

                    system.part.clear()

    def test_d_vel_diff(self):

        system = self.system
        system.part.clear()

        system.lees_edwards.set_params(
            type="steady_shear", velocity=1.5, sheardir=0, shearplanenormal=2)
        system.part.add(id=0, pos=[1, 1, 0.9*system.box_l[2]])
        system.part.add(id=1, pos=[1, 1, 0.1*system.box_l[2]])

        vel_diff = system.velocity_difference(system.part[0], system.part[1])
        np.testing.assert_array_equal(vel_diff, [1.5, 0, 0])

        system.part.clear()
        system.lees_edwards.set_params(type='off')

        system.part.add(id=0, pos=[1, 1, 0.5], v=[0.4, 0.4, 0.4])
        system.part.add(id=1, pos=[1, 1, 1.0], v=[0.1, 0.1, 0.1])
        vel_diff = system.velocity_difference(system.part[0], system.part[1])
        np.testing.assert_array_equal(vel_diff, system.part[0].v-system.part[1].v)

        system.part.clear()

    def test_e_VerletLists(self):
        """ It is not possible to use Verlet lists with Lees Edwards in the current implementation
        we are checking if the correct runtime error message is raised"""

        system = self.system
        system.cell_system.set_domain_decomposition(use_verlet_lists = True)
        system.lees_edwards.set_params(
            type="steady_shear", velocity=1.2, sheardir=2, shearplanenormal=0)

        with self.assertRaises(Exception):
            system.integrator.run(steps=10)

if __name__ == "__main__":
    ut.main()
