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
import espressomd.bond_breakage
import espressomd.interactions
import espressomd.lees_edwards
import espressomd.virtual_sites

import unittest as ut
import unittest_decorators as utx
import tests_common
import numpy as np
import itertools


np.random.seed(42)
params_lin = {'initial_pos_offset': 0.1, 'time_0': 0.1, 'shear_velocity': 1.2}
params_osc = {'initial_pos_offset': 0.1, 'time_0': -2.1, 'amplitude': 2.3,
              'omega': 2.51}
lin_protocol = espressomd.lees_edwards.LinearShear(**params_lin)


def get_lin_pos_offset(time, initial_pos_offset=None,
                       time_0=None, shear_velocity=None):
    return initial_pos_offset + (time - time_0) * shear_velocity


osc_protocol = espressomd.lees_edwards.OscillatoryShear(**params_osc)
off_protocol = espressomd.lees_edwards.Off()


def axis(coord):
    """Return Cartesian axis from coordinate index"""
    res = np.zeros(3)
    mapping = {"x": 0, "y": 1, "z": 2}
    res[mapping[coord]] = 1
    return res


class LeesEdwards(ut.TestCase):

    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists=True)

    time_step = 0.5
    system.time_step = time_step
    direction_permutations = list(itertools.permutations(["x", "y", "z"], 2))

    def setUp(self):
        self.system.time = 0.0

    def tearDown(self):
        system = self.system
        system.part.clear()
        system.lees_edwards.protocol = None
        if espressomd.has_features("VIRTUAL_SITES"):
            system.virtual_sites = espressomd.virtual_sites.VirtualSitesOff()
        if espressomd.has_features("COLLISION_DETECTION"):
            system.collision_detection.set_params(mode="off")

    def test_00_is_none_by_default(self):

        system = self.system

        # Protocol should be off by default
        self.assertIsNone(system.lees_edwards.protocol)
        self.assertIsNone(system.lees_edwards.shear_direction)
        self.assertIsNone(system.lees_edwards.shear_plane_normal)
        self.assertEqual(system.lees_edwards.shear_velocity, 0.)
        self.assertEqual(system.lees_edwards.pos_offset, 0.)

    def test_protocols(self):
        """Test shear velocity and pos offset vs protocol and time"""

        # Linear shear
        system = self.system
        system.time = 3.15
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=lin_protocol)

        # check protocol assignment
        self.assertEqual(system.lees_edwards.protocol, lin_protocol)

        # check parameter setter/getter consistency
        self.assertEqual(
            system.lees_edwards.protocol.get_params(), params_lin)

        # Check pos offset and shear velocity
        expected_pos = params_lin['initial_pos_offset'] + \
            params_lin['shear_velocity'] * (system.time - params_lin['time_0'])
        expected_vel = params_lin['shear_velocity']
        self.assertAlmostEqual(
            system.lees_edwards.shear_velocity, expected_vel)
        self.assertAlmostEqual(
            system.lees_edwards.pos_offset, expected_pos)

        # Check if the offset is determined correctly

        system.time = 0.0

        # Oscillatory shear
        system.lees_edwards.protocol = osc_protocol

        # check parameter setter/getter consistency
        self.assertEqual(system.lees_edwards.protocol.get_params(), params_osc)

        # check pos offset and shear velocity at different times,
        # check that LE offsets are recalculated on simulation time change
        for time in [0., 2.3]:
            system.time = time
            expected_pos = params_osc['initial_pos_offset'] + params_osc['amplitude'] * \
                np.sin(params_osc['omega'] * (time - params_osc['time_0']))
            expected_vel = params_osc['amplitude'] * params_osc['omega'] * \
                np.cos(params_osc['omega'] * (time - params_osc['time_0']))
            self.assertAlmostEqual(
                system.lees_edwards.shear_velocity, expected_vel)
            self.assertAlmostEqual(
                system.lees_edwards.pos_offset, expected_pos)

        # Check that time change during integration updates offsets
        system.integrator.run(1)
        time = system.time
        expected_pos = params_osc['initial_pos_offset'] + params_osc['amplitude'] * \
            np.sin(params_osc['omega'] * (time - params_osc['time_0']))
        expected_vel = params_osc['amplitude'] * params_osc['omega'] * \
            np.cos(params_osc['omega'] * (time - params_osc['time_0']))
        self.assertAlmostEqual(
            system.lees_edwards.shear_velocity,
            expected_vel)
        self.assertAlmostEqual(
            system.lees_edwards.pos_offset, expected_pos)

        # Check that LE shear axes can be overriden
        self.assertEqual(system.lees_edwards.shear_direction, "x")
        self.assertEqual(system.lees_edwards.shear_plane_normal, "y")
        system.lees_edwards.set_boundary_conditions(
            shear_direction="y", shear_plane_normal="z", protocol=lin_protocol)
        self.assertEqual(system.lees_edwards.shear_direction, "y")
        self.assertEqual(system.lees_edwards.shear_plane_normal, "z")

        # Check that LE is disabled correctly via parameter
        system.lees_edwards.protocol = None
        self.assertIsNone(system.lees_edwards.shear_direction)
        self.assertIsNone(system.lees_edwards.shear_plane_normal)
        self.assertEqual(system.lees_edwards.shear_velocity, 0.)
        self.assertEqual(system.lees_edwards.pos_offset, 0.)

        # Check that LE is disabled correctly via boundary conditions setter
        system.lees_edwards.set_boundary_conditions(
            shear_direction="z", shear_plane_normal="y", protocol=lin_protocol)
        system.lees_edwards.set_boundary_conditions(
            shear_direction=99, shear_plane_normal=99, protocol=None)
        self.assertIsNone(system.lees_edwards.shear_direction)
        self.assertIsNone(system.lees_edwards.shear_plane_normal)
        self.assertEqual(system.lees_edwards.shear_velocity, 0.)
        self.assertEqual(system.lees_edwards.pos_offset, 0.)

        # Check that when LE is disabled, protocols can only be
        # initialized via boundary conditions setter, because the
        # shear direction and shear normal must be known
        with self.assertRaisesRegex(RuntimeError, "must be initialized together with 'protocol' on first activation"):
            system.lees_edwards.protocol = lin_protocol

        # Check assertions
        for invalid in ("xy", "t"):
            with self.assertRaisesRegex(ValueError, "Parameter 'shear_direction' is invalid"):
                system.lees_edwards.set_boundary_conditions(
                    shear_direction=invalid, shear_plane_normal="x",
                    protocol=lin_protocol)
            with self.assertRaisesRegex(ValueError, "Parameter 'shear_plane_normal' is invalid"):
                system.lees_edwards.set_boundary_conditions(
                    shear_direction="x", shear_plane_normal=invalid,
                    protocol=lin_protocol)
        for valid in "xyz":
            with self.assertRaisesRegex(ValueError, "Parameters 'shear_direction' and 'shear_plane_normal' must differ"):
                system.lees_edwards.set_boundary_conditions(
                    shear_direction=valid, shear_plane_normal=valid,
                    protocol=lin_protocol)

    def test_boundary_crossing_lin(self):
        """
        A particle crosses the upper and lower boundary with linear shear.
        Check that position and velocity are updated correctly.
        """

        system = self.system

        box_l = np.copy(system.box_l)
        tol = 1e-10

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.set_boundary_conditions(
                shear_direction=shear_direction,
                shear_plane_normal=shear_plane_normal, protocol=lin_protocol)

            system.time = 0.0
            shear_axis = axis(shear_direction)
            pos = box_l - 0.01
            vel = np.array([1., 1., 1.])
            p1 = system.part.add(pos=pos)
            p2 = system.part.add(pos=pos + box_l)

            # Test crossing the upper boundary
            p1.v = vel
            p2.v = vel
            expected_pos = (pos + vel * system.time_step -
                            system.lees_edwards.pos_offset * shear_axis)
            expected_vel = vel - params_lin['shear_velocity'] * shear_axis

            system.integrator.run(1)

            np.testing.assert_allclose(np.copy(p1.v), expected_vel, atol=tol)
            np.testing.assert_allclose(np.copy(p2.v), expected_vel, atol=tol)
            np.testing.assert_allclose(np.copy(p1.pos), expected_pos, atol=tol)
            np.testing.assert_allclose(
                np.copy(p2.pos), expected_pos + box_l, atol=tol)
            np.testing.assert_allclose(
                np.copy(p1.pos_folded), expected_pos - box_l, atol=tol)
            np.testing.assert_allclose(
                np.copy(p2.pos_folded), expected_pos - box_l, atol=tol)

            # Test crossing the lower boundary
            p1.v = -vel
            p2.v = -vel
            expected_pos = p1.pos + p1.v * system.time_step + \
                system.lees_edwards.pos_offset * shear_axis
            system.integrator.run(1)
            np.testing.assert_allclose(np.copy(p1.v), -expected_vel, atol=tol)
            np.testing.assert_allclose(np.copy(p2.v), -expected_vel, atol=tol)
            np.testing.assert_allclose(np.copy(p1.pos), expected_pos, atol=tol)
            np.testing.assert_allclose(
                np.copy(p2.pos), expected_pos + box_l, atol=tol)
            np.testing.assert_allclose(
                np.copy(p1.pos_folded), np.mod(expected_pos - box_l, box_l),
                atol=tol)
            np.testing.assert_allclose(
                np.copy(p2.pos_folded), np.mod(expected_pos - box_l, box_l),
                atol=tol)

    def test_boundary_crossing_off(self):
        """
        A particle crosses the upper and lower boundary without shear.
        Check that position and velocity are unaffected.
        """

        system = self.system
        box_l = np.copy(system.box_l)
        tol = 1e-10

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.set_boundary_conditions(
                shear_direction=shear_direction,
                shear_plane_normal=shear_plane_normal, protocol=off_protocol)

            system.time = 0.0
            pos = box_l - 0.01
            vel = np.array([1., 1., 1.])
            p = system.part.add(pos=pos)

            # Test crossing the upper boundary
            p.v = vel
            expected_pos = pos + vel * system.time_step
            expected_vel = vel
            system.integrator.run(1)
            np.testing.assert_allclose(np.copy(p.v), expected_vel, atol=tol)
            np.testing.assert_allclose(np.copy(p.pos), expected_pos, atol=tol)

            # Test crossing the lower boundary
            p.v = -vel
            expected_pos = p.pos + p.v * system.time_step
            system.integrator.run(1)
            np.testing.assert_allclose(np.copy(p.v), -expected_vel, atol=tol)
            np.testing.assert_allclose(np.copy(p.pos), expected_pos, atol=tol)

    def test_trajectory_reconstruction(self):
        system = self.system
        system.time = 3.4

        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=lin_protocol)

        pos = system.box_l - 0.01
        vel = np.array([0, 1, 0])
        p = system.part.add(pos=pos, v=vel)

        crossing_time = system.time
        system.integrator.run(1)
        np.testing.assert_almost_equal(
            p.lees_edwards_offset, 
            get_lin_pos_offset(crossing_time, **params_lin))
        np.testing.assert_almost_equal(p.lees_edwards_flag, -1)

        system.integrator.run(1)  # no boundary crossing
        np.testing.assert_almost_equal(
            p.lees_edwards_offset, 
            get_lin_pos_offset(crossing_time, **params_lin))

        np.testing.assert_almost_equal(p.lees_edwards_flag, 0)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_distance_vel_diff(self):
        """
        Check distance and velocity difference calculation across LE boundary.
        """

        system = self.system
        epsilon = 0.01

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.set_boundary_conditions(
                shear_direction=shear_direction,
                shear_plane_normal=shear_plane_normal, protocol=lin_protocol)

            shear_axis = axis(shear_direction)
            p1 = system.part.add(
                pos=[epsilon] * 3, v=np.random.random(3), fix=[True] * 3)
            p2 = system.part.add(
                pos=system.box_l - epsilon, v=np.random.random(3), fix=[True] * 3)
            r_euclid = -2 * np.array([epsilon] * 3)

            # check distance
            np.testing.assert_allclose(
                np.copy(system.distance_vec(p1, p2)),
                r_euclid + system.lees_edwards.pos_offset * shear_axis)
            np.testing.assert_allclose(
                np.copy(system.distance_vec(p1, p2)),
                -np.copy(system.distance_vec(p2, p1)))

            # Check velocity difference
            np.testing.assert_allclose(
                np.copy(system.velocity_difference(p1, p2)),
                np.copy(p2.v - p1.v) - system.lees_edwards.shear_velocity * shear_axis)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES", "SOFT_SPHERE"])
    def test_interactions(self):
        """
        We place two particles crossing a boundary and connect them with an
        unbonded and bonded interaction. We test with the resulting stress
        tensor if the offset is included properly.
        """

        system = self.system
        epsilon = 0.01

        for shear_direction, shear_plane_normal in self.direction_permutations:
            system.lees_edwards.set_boundary_conditions(
                shear_direction=shear_direction,
                shear_plane_normal=shear_plane_normal, protocol=lin_protocol)

            p1 = system.part.add(pos=[epsilon] * 3, fix=[True] * 3)
            p2 = system.part.add(pos=system.box_l - epsilon, fix=[True] * 3)

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
                np.copy(system.analysis.pressure_tensor()["bonded"]),
                np.outer(r_12, np.copy(p2.f)) / system.volume())

            np.testing.assert_almost_equal(
                system.analysis.energy()["bonded"],
                0.5 * k_bond * np.copy(system.distance(p1, p2))**2)
            p1.bonds = []

            # Check non-bonded interaction
            k_non_bonded = 3.2
            # NOTE: The force is k*n *distance, hence the 1/2
            system.non_bonded_inter[0, 0].soft_sphere.set_params(
                a=k_non_bonded / 2, n=-2, cutoff=r_cut)
            system.integrator.run(0)
            r_12 = system.distance_vec(p1, p2)

            np.testing.assert_allclose(
                k_non_bonded * r_12, np.copy(p1.f))
            np.testing.assert_allclose(np.copy(p1.f), -np.copy(p2.f))

            np.testing.assert_allclose(
                np.copy(system.analysis.pressure_tensor()["non_bonded"]),
                np.outer(r_12, p2.f) / system.volume())

            np.testing.assert_almost_equal(
                system.analysis.energy()["non_bonded"],
                0.5 * k_non_bonded * np.copy(system.distance(p1, p2))**2)

            system.non_bonded_inter[0, 0].soft_sphere.set_params(
                a=0, n=-2, cutoff=r_cut)
            system.part.clear()

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES", "VIRTUAL_SITES_RELATIVE"])
    def test_virt_sites(self):
        """
        Test placement and force transfer for virtual sites across LE
        boundaries.
        """
        system = self.system
        system.min_global_cut = 2.5
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        tol = 1e-10

        # Construct pair of VS across normal boundary
        system.lees_edwards.protocol = None
        p1 = system.part.add(pos=(2.5, 0.0, 2.5), rotation=[
                             False] * 3, id=0, v=np.array((-1, 2, 3)))
        p2 = system.part.add(pos=(2.5, 1.0, 2.5))
        p2.vs_auto_relate_to(p1)
        p3 = system.part.add(pos=(2.5, 4.0, 2.5))
        p3.vs_auto_relate_to(p1)

        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=lin_protocol)
        # Test position and velocity of VS with Le shift
        old_p3_pos = p3.pos
        expected_p3_pos = old_p3_pos - \
            np.array((get_lin_pos_offset(system.time, **params_lin), 0, 0))
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(np.copy(p3.pos_folded), expected_p3_pos)
        np.testing.assert_allclose(
            p3.v, p1.v + np.array((params_lin["shear_velocity"], 0, 0)))

        # Check distances
        np.testing.assert_allclose(
            np.copy(system.distance_vec(p3, p2)), [0, 2, 0], atol=tol)
        np.testing.assert_allclose(
            np.copy(system.distance_vec(p2, p3)), [0, -2, 0], atol=tol)
        np.testing.assert_allclose(
            np.copy(system.velocity_difference(p3, p2)), [0, 0, 0], atol=tol)
        np.testing.assert_allclose(
            np.copy(system.velocity_difference(p2, p3)), [0, 0, 0], atol=tol)
        system.integrator.run(0)
        np.testing.assert_allclose(
            np.copy(system.distance_vec(p3, p2)), [0, 2, 0], atol=tol)
        system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(system.distance_vec(p3, p2)), [0, 2, 0], atol=tol)

        # Check that force back-transfer matches distance
        p2.ext_force = [1, 0, 0]
        p3.ext_force = -p2.ext_force
        system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p1.torque_lab), [0, 0, -2], atol=tol)

    @utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE", "ROTATION", "DPD"])
    def test_virt_sites_interaction(self):
        """
        A virtual site interacts with a real particle via a DPD interaction
        to get a velocity-dependent force. First we measure a force within the
        primary simulation box as reference. Then we compare it first with the
        situation of the interaction across the Lees Edwards boundary and
        second with the vector from the real particle to the virtual site
        crossing the boundary.
        """

        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

        system.thermostat.set_dpd(kT=0.0, seed=1)
        system.non_bonded_inter[11, 11].dpd.set_params(
            weight_function=0, gamma=1.75, r_cut=2.,
            trans_weight_function=0, trans_gamma=1.5, trans_r_cut=2.0)

        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=2.0, initial_pos_offset=0.0)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)
        p1 = system.part.add(pos=[2.5, 2.5, 2.5], type=10,
                             rotation=3 * (True,), v=(0.0, -0.1, -0.25))
        p2 = system.part.add(pos=(2.5, 3.5, 2.5), type=11)
        p2.vs_auto_relate_to(p1)

        p3 = system.part.add(pos=(2.5, 4.5, 2.5), type=11, v=(2.0, 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        f_p1 = np.copy(p1.f)
        f_p2 = np.copy(p2.f)
        f_p3 = np.copy(p3.f)

        system.part.clear()

        p1 = system.part.add(
            pos=[2.5, 3.75, 2.5], type=10, v=(0.0, -0.1, -0.25),
            rotation=[True] * 3)
        p2 = system.part.add(pos=(2.5, 4.75, 2.5), type=11)
        p2.vs_auto_relate_to(p1)

        p3 = system.part.add(pos=(2.5, 5.75, 2.5), type=11, v=(0.0, 1., 1.25))

        system.integrator.run(0, recalc_forces=True)

        np.testing.assert_array_almost_equal(np.copy(p1.f), f_p1)
        np.testing.assert_array_almost_equal(np.copy(p2.f), f_p2)
        np.testing.assert_array_almost_equal(np.copy(p3.f), f_p3)

        system.part.clear()

        p1 = system.part.add(pos=(2.5, 4.5, 2.5), type=10, v=(0.0, -0.1, -0.25),
                             rotation=3 * (True,))
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

    @utx.skipIfMissingFeatures(
        ["EXTERNAL_FORCES", "VIRTUAL_SITES_RELATIVE", "COLLISION_DETECTION"])
    def test_le_colldet(self):
        system = self.system
        system.min_global_cut = 1.2
        system.time = 0
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=-1.0, initial_pos_offset=0.0)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)

        col_part1 = system.part.add(
            pos=(2.5, 4.5, 2.5), type=30, fix=[True, True, True])
        col_part2 = system.part.add(
            pos=(1.5, 0.5, 2.5), type=30, fix=[True, True, True])

        harm = espressomd.interactions.HarmonicBond(k=1.0, r_0=0.0)
        system.bonded_inter.add(harm)
        virt = espressomd.interactions.Virtual()
        system.bonded_inter.add(virt)

        system.collision_detection.set_params(
            mode="bind_centers", distance=1., bond_centers=harm)

        # After two integration steps we should not have a bond,
        # as the collision detection uses the distant calculation
        # of the short range loop
        system.integrator.run(2)
        bond_list = col_part1.bonds + col_part2.bonds
        np.testing.assert_array_equal(len(bond_list), 0)

        # Bond should be formed on the third integration step
        system.integrator.run(1)

        # One particle should have the bond now.
        bond_list = col_part1.bonds + col_part2.bonds
        np.testing.assert_array_equal(len(bond_list), 1)

        system.part.clear()
        system.collision_detection.set_params(mode="off")

        system.time = 0
        system.lees_edwards.protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=-1.0, initial_pos_offset=0.0)

        system.collision_detection.set_params(
            mode="bind_at_point_of_collision", distance=1., bond_centers=virt,
            bond_vs=harm, part_type_vs=31, vs_placement=1 / 3)

        col_part1 = system.part.add(
            pos=(2.5, 4.5, 2.5), type=30, fix=[True, True, True])
        col_part2 = system.part.add(
            pos=(1.5, 0.5, 2.5), type=30, fix=[True, True, True])

        system.integrator.run(2)
        # We need the distance vector to calculate the positions of the
        # generated VS

        # No bonds present
        bond_list = []
        for p in system.part:
            bond_list = np.append(p.bonds, bond_list)
        np.testing.assert_array_equal(len(bond_list), 0)

        system.integrator.run(1)

        # After the collision detection, we should have four particles
        np.testing.assert_array_equal(len(system.part), 4)

        # Two bonds present
        bond_list = []
        for p in system.part:
            bond_list += p.bonds
        np.testing.assert_array_equal(len(bond_list), 2)

        # We can check on the positions of the generated VS with using the
        # distance of the particles at the time of collision (mi_vec = 1.0).
        # With the placements parameter 1/3 we know the y-component of the
        # generated VS. The other components are inherited from the real
        # particles.
        box_l = np.copy(system.box_l)
        p_vs = system.part.select(virtual=True)
        np.testing.assert_array_almost_equal(
            np.minimum(np.abs(p_vs.pos[:, 0] - col_part1.pos[0]),
                       np.abs(p_vs.pos[:, 0] - col_part2.pos[0])), 0.)
        np.testing.assert_array_almost_equal(p_vs.pos[:, 2], col_part1.pos[2])
        np.testing.assert_array_almost_equal(
            np.sort(np.fmod(p_vs.pos[:, 1] + box_l[1], box_l[1])),
            [1. / 6., box_l[1] - 1. / 6.])

    @utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE"])
    def test_le_breaking_bonds(self):
        system = self.system
        system.min_global_cut = 1.2
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=-1.0, initial_pos_offset=0.0)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)

        harm = espressomd.interactions.HarmonicBond(
            k=1.0, r_0=0.0, r_cut=np.sqrt(2.))
        system.bonded_inter.add(harm)

        p1 = system.part.add(pos=(2.5, 4.5, 2.5))
        p2 = system.part.add(pos=(2.5, 0.5, 2.5))
        p1.add_bond((harm, p2))

        system.bond_breakage[harm] = espressomd.bond_breakage.BreakageSpec(
            breakage_length=np.sqrt(2.), action_type="delete_bond")

        system.integrator.run(3)

        # Bond list should be empty
        bond_list = []
        for p in system.part:
            bond_list += p.bonds
        np.testing.assert_array_equal(len(bond_list), 0)

        system.part.clear()

        # Setup of particles and vs to check if the bind at point of collision
        # will be reverted properly

        harm2 = espressomd.interactions.HarmonicBond(
            k=1.0, r_0=0.0, r_cut=np.sqrt(2.) / 2.)
        system.bonded_inter.add(harm2)

        p1 = system.part.add(pos=(2.5, 4.5, 2.5))
        p2 = system.part.add(pos=(2.5, 0.5, 2.5))

        p3 = system.part.add(pos=(2.5, 4.75, 2.5))
        p3.vs_auto_relate_to(p1)
        p4 = system.part.add(pos=(2.5, 0.25, 2.5))
        p4.vs_auto_relate_to(p2)

        p3.add_bond((harm2, p4))

        system.bond_breakage[harm2] = espressomd.bond_breakage.BreakageSpec(
            breakage_length=np.sqrt(2.) / 2.,
            action_type="revert_bind_at_point_of_collision")

        system.integrator.run(3)

        # Check that all bonds have been removed from the system
        # So the bond list should be empty
        bond_list = []
        for p in system.part:
            bond_list += p.bonds
        np.testing.assert_array_equal(len(bond_list), 0)

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

        np.random.seed(42)
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

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_zz_lj(self):
        """
        Simulate an LJ liquid under linear shear and verify forces.
        This is to make sure that no pairs get lost or are outdated
        in the short range loop. To have deterministic forces, velocity
        capping is used rather than a thermostat.
        """
        system = self.system
        self.setup_lj_liquid()
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=0.3, initial_pos_offset=0.01)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="z", shear_plane_normal="x", protocol=protocol)
        system.integrator.run(1, recalc_forces=True)
        tests_common.check_non_bonded_loop_trace(self, system)

        # Rewind the clock to get back the LE offset applied during force calc
        system.time = system.time - system.time_step
        tests_common.verify_lj_forces(system, 1E-7)

        system.thermostat.set_langevin(kT=.1, gamma=5, seed=2)
        system.integrator.run(50)
        tests_common.check_non_bonded_loop_trace(self, system)


if __name__ == "__main__":
    ut.main()
