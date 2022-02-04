#
# Copyright (C) 2021-2022 The ESPResSo project
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
import espressomd
import espressomd.interactions
import espressomd.lees_edwards as lees_edwards
from espressomd.virtual_sites import VirtualSitesRelative, VirtualSitesOff
from tests_common import verify_lj_forces, check_non_bonded_loop_trace

import unittest as ut
import numpy as np
import itertools


params_lin = {'shear_velocity': 1.2, 'initial_pos_offset': 0.1, 'time_0': 0.1}
lin_protocol = lees_edwards.LinearShear(**params_lin)
params_osc = {'amplitude': 2.3, 'omega': 2.51, 'time_0': -2.1}


def axis(coord):
    """Return Cartesian axis from coordinate index"""
    res = np.zeros(3)
    res[coord] = 1
    return res


class LeesEdwards(ut.TestCase):

    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists=True)

    time_step = 0.5
    system.time_step = time_step
    direction_permutations = list(itertools.permutations([0, 1, 2], 2))

    def setUp(self):
        self.system.time = 0.0

    def tearDown(self):
        self.system.part.clear()
        self.system.lees_edwards.protocol = None

    def test_00_is_off_by_default(self):

        system = self.system

        # Protocol should be off by default
        self.assertIsNone(system.lees_edwards.protocol)
        self.assertEqual(system.lees_edwards.shear_velocity, 0)
        self.assertEqual(system.lees_edwards.pos_offset, 0)

    def test_protocols(self):
        """Tests shear velocity and pos offset vs protocol and time"""

        # Linear shear
        system = self.system
        system.time = 3.15
        system.lees_edwards.protocol = lin_protocol

        # check protocol assignment
        self.assertEqual(self.system.lees_edwards.protocol, lin_protocol)

        # check parameter setter/getter consistency
        self.assertEqual(
            self.system.lees_edwards.protocol.get_params(),
            params_lin)

        # Check pos offset and shear velocity
        expected_pos = params_lin['initial_pos_offset'] + \
            params_lin['shear_velocity'] * (system.time - params_lin['time_0'])
        expected_vel = params_lin['shear_velocity']
        self.assertAlmostEqual(
            system.lees_edwards.shear_velocity,
            expected_vel)
        self.assertAlmostEqual(
            system.lees_edwards.pos_offset, expected_pos)

        # Check if the offset is determined correctly

        system.time = 0.0

        # Oscillatory shear
        system.lees_edwards.protocol = lees_edwards.OscillatoryShear(
            **params_osc)

        # check parameter setter/getter consistency
        self.assertEqual(
            self.system.lees_edwards.protocol.get_params(),
            params_osc)

        # check pos offset and shear velocity at different times,
        # check that LE offsets are recalculated on simulation tmie change
        for time in [0., 2.3]:
            system.time = time
            expected_pos = params_osc['amplitude'] * \
                np.sin(params_osc['omega'] * (time - params_osc['time_0']))
            expected_vel = params_osc['amplitude'] * params_osc['omega'] * \
                np.cos(params_osc['omega'] * (time - params_osc['time_0']))
            self.assertAlmostEqual(
                system.lees_edwards.shear_velocity,
                expected_vel)
            self.assertAlmostEqual(
                system.lees_edwards.pos_offset, expected_pos)
        # Check that time change during integration updates offsets
        system.integrator.run(1)
        time = system.time
        expected_pos = params_osc['amplitude'] * \
            np.sin(params_osc['omega'] * (time - params_osc['time_0']))
        expected_vel = params_osc['amplitude'] * params_osc['omega'] * \
            np.cos(params_osc['omega'] * (time - params_osc['time_0']))
        self.assertAlmostEqual(
            system.lees_edwards.shear_velocity,
            expected_vel)
        self.assertAlmostEqual(
            system.lees_edwards.pos_offset, expected_pos)

        # Check that LE is diabled correctly
        system.lees_edwards.protocol = None
        self.assertEqual(
            system.lees_edwards.shear_velocity, 0)
        self.assertEqual(
            system.lees_edwards.pos_offset, 0)

    def test_boundary_crossing(self):
        """A particle crosses the upper and lower boundary to test if position
        and velocity are updated correctly."""

        system = self.system

        # Set up a one particle system and check the position offset after crossing the boundary
        # Test for upper boundary

        system.lees_edwards.protocol = lin_protocol

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.shear_direction = shear_direction
            system.lees_edwards.shear_plane_normal = shear_plane_normal

            shear_axis = axis(shear_direction)

            pos = system.box_l - 0.01
            vel = np.array([1, 1, 1])
            p = system.part.add(pos=pos, v=vel)

            system.time = 0.0
            expected_pos = (pos + vel * system.time_step -
                            system.lees_edwards.pos_offset * shear_axis)
            expected_delta_vel = -params_lin['shear_velocity'] * shear_axis

            system.integrator.run(1)

            np.testing.assert_almost_equal(p.v, vel + expected_delta_vel)
            np.testing.assert_almost_equal(p.pos, expected_pos)

            # Test the crossing lower boundary
            p.v = -vel
            expected_pos = p.pos + p.v * system.time_step + \
                system.lees_edwards.pos_offset * shear_axis
            system.integrator.run(1)
            np.testing.assert_almost_equal(p.pos, expected_pos)
            np.testing.assert_almost_equal(p.v, -vel - expected_delta_vel)

    def test_trajectory_reconstruction(self):
        system = self.system

        system.lees_edwards.protocol = lees_edwards.LinearShear(
            shear_velocity=1., initial_pos_offset=0.0, time_0=0.0)
        system.lees_edwards.shear_direction = 0
        system.lees_edwards.shear_plane_normal = 1

        pos = system.box_l - 0.01
        vel = np.array([0, 1, 0])
        p = system.part.add(pos=pos, v=vel)

        system.integrator.run(1)

        np.testing.assert_almost_equal(
            p.lees_edwards_flag * 1.0 * system.time_step * 0.5,
            p.lees_edwards_offset)
        np.testing.assert_almost_equal(-1, p.lees_edwards_flag)

        offset1 = p.lees_edwards_flag * 1.0 * system.time_step * 0.5

        system.integrator.run(1)

        np.testing.assert_almost_equal(
            offset1 - 1.0 * 0.5, p.lees_edwards_offset)
        np.testing.assert_almost_equal(0, p.lees_edwards_flag)


#        x1 = pos[0]
#        y1 = pos[1]
#
#        x2 = p.pos[0]
#        y2 = p.pos[1]
#
#        x3 = pos[0] + p.lees_edwards_flag * p.lees_edwards_offset
#        y3 = system.box_l[1]
#
#        slope1 = (y2-y1)/(x2-x1)
#        print("Slope:", slope1)
#
#        slope2 = (y3-y1)/(x3-x1)
#        print("Slope:", slope2)
#
#        plt.plot([x1,x2,x3],[y1,y2,y3], 's')
#        plt.show()


    def test_distance_vel_diff(self):
        """check distance and velocity difference calculation across LE boundary
        """

        system = self.system
        system.lees_edwards.protocol = lin_protocol
        epsilon = 0.01

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.shear_direction = shear_direction
            system.lees_edwards.shear_plane_normal = shear_plane_normal
            shear_axis = axis(shear_direction)

            system.lees_edwards.protocol = lin_protocol
            p1 = system.part.add(
                pos=[epsilon] * 3, v=np.random.random(3), fix=[1] * 3)
            p2 = system.part.add(
                pos=system.box_l - epsilon, v=np.random.random(3), fix=[1] * 3)
            r_euclid = -2 * np.array([epsilon] * 3)

            # check distance
            np.testing.assert_allclose(
                system.distance_vec(p1, p2), r_euclid + system.lees_edwards.pos_offset * shear_axis)
            np.testing.assert_allclose(
                np.copy(system.distance_vec(p1, p2)),
                -np.copy(system.distance_vec(p2, p1)))

            # Check velocity difference
            np.testing.assert_allclose(
                system.velocity_difference(p1, p2),
                p2.v - p1.v - system.lees_edwards.shear_velocity * shear_axis)

    def test_interactions(self):
        """We place two particles crossing a boundary and connect them with an unbonded
           and bonded interaction. We test with the resulting stress tensor if the offset
           is included properly"""

        system = self.system
        system.lees_edwards.protocol = lin_protocol
        epsilon = 0.01

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.shear_direction = shear_direction
            system.lees_edwards.shear_plane_normal = shear_plane_normal

            system.lees_edwards.protocol = lin_protocol
            p1 = system.part.add(pos=[epsilon] * 3, fix=[1] * 3)
            p2 = system.part.add(pos=system.box_l - epsilon, fix=[1] * 3)

            # check bonded interaction
            k_bond = 2.1
            r_cut = 1.5
            harm = espressomd.interactions.HarmonicBond(k=k_bond, r_0=0.0)
            system.bonded_inter.add(harm)
            p1.add_bond((harm, p2))

            r_12 = np.copy(system.distance_vec(p1, p2))
            system.integrator.run(0)

            np.testing.assert_allclose(np.copy(p1.f), k_bond * r_12)
            np.testing.assert_allclose(np.copy(p1.f), -np.copy(p2.f))

            np.testing.assert_allclose(
                system.analysis.pressure_tensor()["bonded"],
                np.outer(r_12, np.copy(p2.f)) / system.volume())

            np.testing.assert_almost_equal(
                system.analysis.energy()["bonded"],
                1 / 2 * k_bond * system.distance(p1, p2)**2)
            p1.bonds = []

            # Check non bonded interaction
            k_non_bonded = 3.2
            # NOTE: The force is k*n *distance, hence the 1/2
            system.non_bonded_inter[0, 0].soft_sphere.set_params(
                a=k_non_bonded / 2, n=-2, cutoff=r_cut)
            system.integrator.run(0)
            r_12 = system.distance_vec(p1, p2)

            np.testing.assert_allclose(
                (k_non_bonded) * r_12,
                np.copy(p1.f))
            np.testing.assert_allclose(np.copy(p1.f), -np.copy(p2.f))

            np.testing.assert_allclose(
                system.analysis.pressure_tensor()["non_bonded"],
                np.outer(r_12, p2.f) / system.volume())

            np.testing.assert_almost_equal(
                system.analysis.energy()["non_bonded"],
                1 / 2 * k_non_bonded * system.distance(p1, p2)**2)

            system.non_bonded_inter[0, 0].soft_sphere.set_params(
                a=0, n=-2, cutoff=r_cut)
            system.part.clear()

    def test_virt_sites(self):
        """Tests placement and force transfer for virtual sites across LE
        boundaries. """
        system = self.system
        system.min_global_cut = 2.5
        self.system.virtual_sites = VirtualSitesRelative()
        tol = 1e-10

        # Construct pair of VS across normal boudnary
        system.lees_edwards.protocol = None
        p1 = system.part.add(id=0, pos=[2.5, 0.0, 2.5], rotation=(0, 0, 0))
        p2 = system.part.add(pos=(2.5, 1.0, 2.5))
        p2.vs_auto_relate_to(p1)
        p3 = system.part.add(pos=(2.5, 4.0, 2.5))
        p3.vs_auto_relate_to(p1)
        system.integrator.run(1)

        system.lees_edwards.protocol = lin_protocol
        system.lees_edwards.shear_direction = 0
        system.lees_edwards.shear_plane_normal = 1
        system.integrator.run(1)
        np.testing.assert_allclose(
            system.distance_vec(p3, p2), [0, 2, 0], atol=tol)
        np.testing.assert_allclose(
            system.velocity_difference(p3, p2), [0, 0, 0], atol=tol)
        system.integrator.run(0)
        np.testing.assert_allclose(
            system.distance_vec(p3, p2), [0, 2, 0], atol=tol)
        system.integrator.run(1)
        np.testing.assert_allclose(
            system.distance_vec(p3, p2), [0, 2, 0], atol=tol)

        # Check that force back transfer matches distance
        p2.ext_force = [1, 0, 0]
        p3.ext_force = -p2.ext_force
        system.integrator.run(1)
        np.testing.assert_allclose(p1.torque_lab, [0, 0, -2], atol=tol)

    def test_virt_sites_interaction(self):
        """A virtual site interacts with a real particle via a DPD interaction to get
           a velocity dependent force. First we measure a force within the primary
           simulation box as reference. Then we compare it first with the situation of
           the interaction across the Lees Edwards boundary and second with the vector
           from the real particle to the virtual site crossing the boundary."""

        system = self.system
        self.system.virtual_sites = VirtualSitesRelative()

        system.thermostat.set_dpd(kT=0.0, seed=1)
        system.non_bonded_inter[11, 11].dpd.set_params(
            weight_function=0, gamma=1.75, r_cut=2.,
            trans_weight_function=0, trans_gamma=1.5, trans_r_cut=2.0)

        system.lees_edwards.protocol = lees_edwards.LinearShear(
            shear_velocity=2.0, initial_pos_offset=0.0)
        system.lees_edwards.shear_direction = 0
        system.lees_edwards.shear_plane_normal = 1
        p1 = system.part.add(
            pos=[2.5, 2.5, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 3.5, 2.5), type=11)
        p2.vs_auto_relate_to(p1)

        p3 = system.part.add(pos=(2.5, 4.5, 2.5), type=11, v=(2.0, 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        f_p1 = np.copy(p1.f)
        f_p2 = np.copy(p2.f)
        f_p3 = np.copy(p3.f)

        system.part.clear()

        p1 = system.part.add(
            pos=[2.5, 3.75, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 4.75, 2.5), type=11)
        p2.vs_auto_relate_to(p1)

        p3 = system.part.add(pos=(2.5, 5.75, 2.5), type=11, v=(0.0, 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
        np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
        np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)

        system.part.clear()

        p1 = system.part.add(
            pos=[2.5, 4.5, 2.5], rotation=(1, 1, 1), type=10, v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 5.5, 2.5), type=11)
        p2.vs_auto_relate_to(p1)

        p3 = system.part.add(pos=(2.5, 6.5, 2.5), type=11, v=(0., 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
        np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
        np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)

        system.part.clear()

        system.non_bonded_inter[11, 11].dpd.set_params(
            weight_function=0, gamma=0, r_cut=0,
            trans_weight_function=0, trans_gamma=0, trans_r_cut=0)
        system.virtual_sites = VirtualSitesOff()

    def setup_lj_liquid(self):
        system = self.system
        system.cell_system.set_n_square(use_verlet_lists=False)
        # Parameters
        n = 100
        phi = 0.4
        sigma = 1.
        eps = 1
        cut = sigma * 2**(1 / 6)

        # box
        l = (n / 6. * np.pi * sigma**3 / phi)**(1. / 3.)

        # Setup
        system.box_l = [l, l, l]
        system.lees_edwards.protocol = None

        system.time_step = 0.01
        system.thermostat.turn_off()

        system.part.add(pos=np.random.random((n, 3)) * l)

        # interactions
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=eps, sigma=sigma, cutoff=cut, shift="auto")
        # Remove overlap
        system.integrator.set_steepest_descent(
            f_max=0, gamma=0.05, max_displacement=0.05)
        while system.analysis.energy()["total"] > 0.5 * n:
            system.integrator.run(5)

        system.integrator.set_vv()

    def test_z5_lj(self):
        """
        Simulates an LJ liquid under linear shear and verifies forces. This is to
        make sure that no pairs get lost or are outdated in the short range loop.
        To have deterministic forces, velocity capping is used rather than a
        thermostat."""
        system = self.system
        self.setup_lj_liquid()
        system.lees_edwards.protocol = lees_edwards.LinearShear(
            shear_velocity=0.3, initial_pos_offset=0.01)
        system.lees_edwards.shear_direction = 2
        system.lees_edwards.shear_plane_normal = 0
        system.integrator.run(1, recalc_forces=True)
        check_non_bonded_loop_trace(system)

        # Rewind the clock to get back the LE offset applied during force calc
        system.time = system.time - system.time_step
        verify_lj_forces(system, 1E-7)

        system.thermostat.set_langevin(kT=.1, gamma=5, seed=2)
        system.integrator.run(50)
        check_non_bonded_loop_trace(system)


if __name__ == "__main__":
    ut.main()
