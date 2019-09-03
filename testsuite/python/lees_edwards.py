from __future__ import print_function

import espressomd
from espressomd.interactions import HarmonicBond
import espressomd.lees_edwards as lees_edwards

import unittest as ut
import numpy as np


@ut.skipIf(not espressomd.has_features(["LEES_EDWARDS"]),
           "Feature not available, skipping test!")
class LeesEdwards(ut.TestCase):

    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists=False)
    system.set_random_state_PRNG()

    time_step = 0.5
    system.time_step = time_step

    def test_a_Protocol(self):
        """This test calculates if the offset and velocity of the Lees Edwards
        function are calculated correctly based on the input variables"""

        system = self.system

        # Protocol should be off by default
        self.assertTrue(self.system.lees_edwards.protocol is None)

        # Check if the set protocol is stored correctly
        params = {
            'shear_velocity': 1.2,
            'shear_direction': 2,
            'shear_plane_normal': 0,
            'initial_pos_offset': 0.1}
        new_protocol = lees_edwards.LinearShear(**params)
        system.lees_edwards.protocol = new_protocol 
        self.assertEqual(self.system.lees_edwards.protocol, new_protocol)
        self.assertEqual(
            self.system.lees_edwards.protocol.get_params(),
            params)

        # Check if the offset is determined correctly

        system.time = 0.0

        frequency = np.sqrt(2.0)
        omega = 2 * np.pi * frequency
        amplitude = 1.3
        sheardir = 0
        shear_plane_normal = 1
        system.lees_edwards.protocol = lees_edwards.OscillatoryShear(
            frequency=omega,
            amplitude=amplitude, 
            shear_direction=sheardir, 
            shear_plane_normal=shear_plane_normal)

        for time in np.arange(0.0, 100.0, system.time_step * 10):

            system.time = time
            offset = amplitude * np.sin(omega * time)
            velocity = omega * amplitude * np.cos(omega * time)
            self.assertAlmostEqual(
                system.lees_edwards.shear_velocity,
                velocity)
            self.assertAlmostEqual(system.lees_edwards.pos_offset, offset)

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
            for shear_plane_normal in dir:
                if sheardir != shear_plane_normal:

                    system.time = 0.0
                    system.lees_edwards.protocol = \
                        lees_edwards.LinearShear(
                            shear_velocity=velocity, 
                          shear_direction=sheardir, 
                          shear_plane_normal=shear_plane_normal)

                    pos = np.full([3], 2.5)
                    pos[shear_plane_normal] = system.box_l[
                        shear_plane_normal] - 0.05
                    vel = np.zeros([3])
                    vel[shear_plane_normal] = 0.2
                    system.part.add(pos=pos, v=vel, id=0, type=0)

                    pos_change = np.zeros([3])
                    pos_change[sheardir] = -system.time_step * 0.5 * velocity
                    pos_change[
                        shear_plane_normal] = vel[
                            shear_plane_normal] * system.time_step
                    vel_change = np.zeros([3])
                    vel_change[sheardir] = -velocity

                    expected_pos = system.part[0].pos + pos_change
                    expected_vel = system.part[0].v + vel_change

                    system.integrator.run(1)

                    np.testing.assert_almost_equal(
                        np.copy(system.part[0].pos), expected_pos)
                    np.testing.assert_almost_equal(
                        np.copy(system.part[0].v), expected_vel)

                    system.part.clear()

                    # Test for lower boundary

                    system.time = 0.0
                    pos = np.full([3], 2.5)
                    pos[shear_plane_normal] = 0.05
                    vel = np.zeros([3])
                    vel[shear_plane_normal] = -0.2
                    system.part.add(pos=pos, v=vel, id=0, type=0)

                    pos_change = np.zeros([3])
                    pos_change[sheardir] = system.time_step * 0.5 * velocity
                    pos_change[
                        shear_plane_normal] = vel[
                            shear_plane_normal] * system.time_step
                    vel_change = np.zeros([3])
                    vel_change[sheardir] = velocity

                    expected_pos = system.part[0].pos + pos_change
                    expected_vel = system.part[0].v + vel_change

                    system.integrator.run(1)

                    np.testing.assert_almost_equal(
                        np.copy(system.part[0].pos), expected_pos)
                    np.testing.assert_almost_equal(
                        np.copy(system.part[0].v), expected_vel)

                    system.part.clear()

    def test_c_interactions(self):
        """We place two particles crossing a boundary and connect them with an unbonded
           and bonded interaction. We test with the resulting stress tensor if the offset
           is included properly"""

        system = self.system
        system.part.clear()

        dir = [0, 1, 2]

        for sheardir in dir:
            for shear_plane_normal in dir:
                if sheardir != shear_plane_normal:

                    offset = 1.0
                    system.lees_edwards.protocol = \
                        lees_edwards.LinearShear(
                            initial_pos_offset=offset,
                            shear_velocity=0,
                            shear_direction=sheardir,
                            shear_plane_normal=shear_plane_normal)

                    pos1 = np.full([3], 2.5)
                    pos1[shear_plane_normal] = 4.75
                    pos2 = np.full([3], 2.5)
                    pos2[shear_plane_normal] = 0.25

                    system.part.add(id=0, pos=pos1, fix=[1, 1, 1])
                    system.part.add(id=1, pos=pos2, fix=[1, 1, 1])

                    r = system.part[1].pos - system.part[0].pos
                    r[sheardir] += offset
                    r[shear_plane_normal] += system.box_l[0]

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
                            analytical_bondedstress[i, j] += k * r[i] * r[j]
                    analytical_bondedstress /= (system.box_l[0]**3.0)

                    analytical_nonbondedstress = np.zeros([3, 3])
                    for i in range(3):
                        for j in range(3):
                            analytical_nonbondedstress[
                                i, j] += 2 * k * r[i] * r[j]
                    analytical_nonbondedstress /= (system.box_l[0]**3.0)

                    np.testing.assert_array_equal(
                        simulated_bondedstress, analytical_bondedstress)
                    np.testing.assert_array_equal(
                        simulated_nonbondedstress, analytical_nonbondedstress)

                    system.part.clear()

    def test_d_vel_diff(self):

        system = self.system
        system.part.clear()

        system.lees_edwards.protocol = \
            lees_edwards.LinearShear(
                shear_velocity=1.5,
                shear_direction=0,
                shear_plane_normal=2)
        p1 = system.part.add(
            id=0, pos=[1, 1, 0.9 * system.box_l[2]], v=[2, -3.2, 3])
        p2 = system.part.add(
            id=1, pos=[1, 1, 0.1 * system.box_l[2]], v=[7, 6.1, -1])

        vel_diff = system.velocity_difference(p1, p2)
        np.testing.assert_array_equal(
            vel_diff, p2.v - p1.v + np.array([system.lees_edwards.shear_velocity, 0, 0]))

        system.part.clear()
        system.lees_edwards.protocol = lees_edwards.Off()

        p1 = system.part.add(id=0, pos=[1, 1, 0.5], v=[0.4, 0.4, 0.4])
        p2 = system.part.add(id=1, pos=[1, 1, 1.0], v=[0.1, 0.1, 0.1])
        vel_diff = system.velocity_difference(p1, p2)
        np.testing.assert_array_equal(
            vel_diff,
           p2.v - p1.v)

        system.part.clear()

    def test_e_VerletLists(self):
        """ It is not possible to use Verlet lists with Lees Edwards in the current implementation
        we are checking if the correct runtime error message is raised"""

        system = self.system
        system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        system.lees_edwards.protocol = \
            lees_edwards.LinearShear(
                shear_velocity=1.2,
                shear_direction=2,
                shear_plane_normal=0)

        with self.assertRaises(Exception):
            system.integrator.run(steps=10)

if __name__ == "__main__":
    ut.main()
