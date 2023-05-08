#
# Copyright (C) 2010-2022 The ESPResSo project
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

import espressomd.lb
import espressomd.shapes
import unittest as ut
import unittest_decorators as utx
import numpy as np
import itertools


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBoundaryVelocityTest(ut.TestCase):
    """
    Various tests to check the interaction of lb velocity boundary conditions and the fluid
    """

    lb_params = {'agrid': 0.6,
                 'density': 0.5,
                 'kinematic_viscosity': 3.2,
                 'tau': 0.7}
    system = espressomd.System(box_l=3 * [8 * lb_params['agrid']])
    system.time_step = lb_params['tau']
    system.cell_system.skin = 0.1

    def tearDown(self):
        self.system.actors.clear()

    def setUp(self):
        self.lb_fluid = espressomd.lb.LBFluidWalberla(**self.lb_params)
        self.system.actors.add(self.lb_fluid)

    def check_wall_slip(self, v_boundary, atol):
        """
        Check that the fluid adopts the velocity set by the boundary conditions.
        """

        agrid = self.lb_params['agrid']
        wall_shape_left = espressomd.shapes.Wall(normal=[1, 0, 0], dist=agrid)
        wall_shape_right = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - agrid))
        for shape in [wall_shape_left, wall_shape_right]:
            self.lb_fluid.add_boundary_from_shape(shape, v_boundary)

        # fluid in contact with moving boundary adopts same velocity
        self.system.integrator.run(200)
        v_fluid = np.copy(self.lb_fluid[2, 1, 3].velocity)
        np.testing.assert_allclose(v_fluid, v_boundary, atol=atol)

        # velocity in the middle needs to propagate first
        self.system.integrator.run(200)
        v_fluid = np.copy(self.lb_fluid[2, 1, 3].velocity)
        np.testing.assert_allclose(v_fluid, v_boundary, atol=atol)

    def test_wall_slip_parallel(self):
        v_boundary = [0, 0, 0.07]
        self.check_wall_slip(v_boundary, 2e-4)

    def test_wall_slip_nonparallel(self):
        v_boundary = [0.03, 0.02, 0.01]
        self.check_wall_slip(v_boundary, 5e-4)

    def test_boundary_readout(self):
        """
        Test the read part of the boundary property of lb nodes.
        """
        v_boundary = [0.03, 0.05, 0.07]
        wall_shape = espressomd.shapes.Wall(
            normal=[1, 0, 0], dist=self.lb_params['agrid'])
        self.lb_fluid.add_boundary_from_shape(wall_shape, v_boundary)

        # check non_boundary node
        bound_cond = self.lb_fluid[4, 4, 4].boundary
        self.assertIsNone(bound_cond)

        # on boundary
        bound_cond = self.lb_fluid[0, 2, 4].boundary
        np.testing.assert_array_almost_equal(bound_cond.velocity, v_boundary)

        # TODO WALBERLA boundary_force

    def test_velocity_bounce_back_class(self):
        """
        Test setters and getters of :class:`espressomd.lb.VelocityBounceBack`
        """
        with self.assertRaises(ValueError):
            bound_cond = espressomd.lb.VelocityBounceBack([1, 2])
        v = [1, 2, 17.4]
        bound_cond = espressomd.lb.VelocityBounceBack(v)
        np.testing.assert_array_almost_equal(bound_cond.velocity, v)

    def test_boundary_setting(self):
        """
        Test setting and un-setting individual lb boundary nodes.
        """
        v_boundary = [0.02, 0.01, 0.03]
        bound_cond = espressomd.lb.VelocityBounceBack(v_boundary)

        with self.assertRaises(TypeError):
            self.lb_fluid[1, 2, 3].boundary = 17
        with self.assertRaises(TypeError):
            self.lb_fluid[1, 2, 3].boundary = np.array([1, 2, 3])

        self.lb_fluid[1, 2, 3].boundary = bound_cond
        np.testing.assert_array_almost_equal(
            self.lb_fluid[1, 2, 3].boundary.velocity, bound_cond.velocity)

        self.lb_fluid[1, 2, 3].boundary = None
        self.assertIsNone(self.lb_fluid[1, 2, 3].boundary)

    def test_nodes_inside_shape_line(self):
        """
        Test if the ``get_nodes_inside_shape`` method correctly identifies
        the grid points inside a line.
        """
        agrid = self.lb_params['agrid']
        cyl = espressomd.shapes.Cylinder(center=agrid * np.array([1.5, 2.5, 5.5]),
                                         axis=[0, 0, 1],
                                         length=2.1 * agrid,
                                         radius=0.5 * agrid)
        nodes_in_boundary = self.lb_fluid.get_nodes_inside_shape(cyl)
        idxs_in_boundary = set(tuple(node.index) for node in nodes_in_boundary)

        idx_ref = {(1, 2, 5), (1, 2, 4), (1, 2, 6)}
        self.assertSetEqual(idxs_in_boundary, idx_ref)

    def test_nodes_inside_shape_cylinder(self):
        """
        Test if the ``get_nodes_inside_shape`` method correctly identifies
        the grid points inside a cylinder.
        """
        agrid = self.lb_params['agrid']
        cyl = espressomd.shapes.Cylinder(center=agrid * np.array([1.5, 1.5, 0.5]),
                                         axis=[0, 0, 1],
                                         length=2.0 * self.system.box_l[2],
                                         radius=2.0 * agrid)
        nodes_in_boundary = self.lb_fluid.get_nodes_inside_shape(cyl)
        idxs_in_boundary = list(tuple(node.index)
                                for node in nodes_in_boundary)

        for node in idxs_in_boundary:
            self.assertIn(node[0], [0, 1, 2])
            self.assertIn(node[1], [0, 1, 2])
            self.assertIn(node[2], np.arange(8))

    def test_nodes_inside_shape_cube(self):
        """
        Test if the ``get_nodes_inside_shape`` method correctly identifies
        the grid points inside a cube.
        """
        agrid = self.lb_params['agrid']
        prism = espressomd.shapes.Rhomboid(a=2 * agrid * np.array([1, 0, 0]),
                                           b=3 * agrid * np.array([0, 1, 0]),
                                           c=4 * agrid * np.array([0, 0, 1]),
                                           corner=agrid * np.array([1, 1, 1]),
                                           direction=1)
        nodes_in_boundary = self.lb_fluid.get_nodes_inside_shape(prism)
        idxs_in_boundary = set(tuple(node.index) for node in nodes_in_boundary)

        idx_ref = set(itertools.product(range(1, 3), range(1, 4), range(1, 5)))
        self.assertSetEqual(idxs_in_boundary, idx_ref)

    def test_shape_bitmask(self):
        """
        Test if the ``get_shape_bitmask`` method correctly identifies the grid
        points inside a shape and matches the LB ``is_boundary`` property.
        """
        def get_masks(shape):
            """
            Get the shape mask and the LB boundary mask.
            """
            self.lb_fluid.add_boundary_from_shape(shape)
            lb_bitmask = np.copy(self.lb_fluid[:, :, :].is_boundary)
            shape_bitmask = self.lb_fluid.get_shape_bitmask(shape)
            self.lb_fluid.clear_boundaries()
            return lb_bitmask.astype(int), shape_bitmask.astype(int)

        agrid = self.lb_params['agrid']

        # check a prism
        for nudge_corner in (0.5 + 1e-6, 1.0, 1.5 - 1e-6):
            shape = espressomd.shapes.Rhomboid(
                a=2 * agrid * np.array([1, 0, 0]),
                b=3 * agrid * np.array([0, 1, 0]),
                c=4 * agrid * np.array([0, 0, 1]),
                corner=agrid * nudge_corner * np.array([1, 1, 1]),
                direction=1)
            lb_bitmask, shape_bitmask = get_masks(shape)
            np.testing.assert_array_equal(shape_bitmask, lb_bitmask)
            np.testing.assert_array_equal(shape_bitmask[1:3, 1:4, 1:5], 1)
            shape_bitmask[1:3, 1:4, 1:5] = 0
            np.testing.assert_array_equal(shape_bitmask, 0)

        # check a sphere
        for nudge_radius in (-0.1, -0.01, 0., 0.01, 0.1):
            for nudge_center in ([0.1, 0., 0.], [0., 0.15, 0.20]):
                shape = espressomd.shapes.Sphere(
                    center=4 * agrid * np.array([1, 1, 1]) + nudge_center,
                    radius=3 * agrid + nudge_radius)
                lb_bitmask, shape_bitmask = get_masks(shape)
                np.testing.assert_array_equal(shape_bitmask, lb_bitmask)

    def test_edge_detection_x(self):
        self.check_edge_detection(0)

    def test_edge_detection_y(self):
        self.check_edge_detection(1)

    def test_edge_detection_z(self):
        self.check_edge_detection(2)

    def check_edge_detection(self, axis):
        """
        Test if the ``edge_detection`` method correctly identifies the grid
        points on the surface of a cube and on the surface of a square
        column (finite or infinite, periodic or aperiodic).
        """
        def get_surface_indices(mask, periodicity):
            idx = espressomd.lb.edge_detection(mask, periodicity)
            return set(map(tuple, idx))

        def roll_product(a, b, c):
            """
            Calculate ``itertools.product`` of 3 objects that are rolled.
            """
            collection = np.array([list(a), list(b), list(c)], dtype=object)
            return itertools.product(*np.roll(collection, axis))

        def create_column_shape_roll(lengths, corner):
            """
            Create a prism with lengths and corner that are rolled.
            """
            lengths = np.roll(lengths, axis)
            corner = np.roll(corner, axis)
            return espressomd.shapes.Rhomboid(
                a=lengths[0] * agrid * np.array([1, 0, 0]),
                b=lengths[1] * agrid * np.array([0, 1, 0]),
                c=lengths[2] * agrid * np.array([0, 0, 1]),
                corner=agrid * corner,
                direction=1)

        agrid = self.lb_params['agrid']
        periodic = np.roll([True, True, True], axis)
        aperiodic = np.roll([False, False, False], axis)

        # check a simple cube
        cube = create_column_shape_roll([4, 4, 4], [1, 1, 1])
        self.lb_fluid.add_boundary_from_shape(cube)
        cube_mask = np.copy(self.lb_fluid[:, :, :].is_boundary.astype(bool))
        idx_ref = set(roll_product(range(1, 5), range(1, 5), range(1, 5)))
        for item in roll_product(range(2, 4), range(2, 4), range(2, 4)):
            idx_ref.remove(item)

        idxs_on_surface = get_surface_indices(cube_mask, periodic)
        self.assertSetEqual(idxs_on_surface, idx_ref)

        self.lb_fluid.clear_boundaries()

        # create an infinite square column
        col = create_column_shape_roll([8, 4, 4], [0, 1, 1])
        self.lb_fluid.add_boundary_from_shape(col)
        col_mask = np.copy(self.lb_fluid[:, :, :].is_boundary.astype(bool))
        idx_ref = set(roll_product(range(0, 8), range(1, 5), range(1, 5)))
        for item in roll_product(range(0, 8), range(2, 4), range(2, 4)):
            idx_ref.remove(item)

        # with periodicity: check neither ends are in contact with fluid
        idxs_on_surface = get_surface_indices(col_mask, periodic)
        self.assertSetEqual(idxs_on_surface, idx_ref)

        # without periodicity: check neither ends are in contact with fluid
        idxs_on_surface = get_surface_indices(col_mask, aperiodic)
        self.assertSetEqual(idxs_on_surface, idx_ref)

        self.lb_fluid.clear_boundaries()

        # create a finite square column; both ends of the columns are in
        # contact with a thin slice of fluid
        col = create_column_shape_roll([7, 4, 4], [0, 1, 1])
        self.lb_fluid.add_boundary_from_shape(col)
        col_mask = np.copy(self.lb_fluid[:, :, :].is_boundary.astype(bool))
        idx_ref = set(roll_product(range(0, 7), range(1, 5), range(1, 5)))

        # with periodicity: check both ends are in contact with fluid
        for item in roll_product(range(1, 6), range(2, 4), range(2, 4)):
            idx_ref.remove(item)
        idxs_on_surface = get_surface_indices(col_mask, periodic)
        self.assertSetEqual(idxs_on_surface, idx_ref)

        # without periodicity: check one end of the column is no longer in
        # contact with the fluid
        for item in roll_product(range(0, 1), range(2, 4), range(2, 4)):
            idx_ref.remove(item)
        idxs_on_surface = get_surface_indices(col_mask, aperiodic)
        self.assertSetEqual(idxs_on_surface, idx_ref)

    def test_calc_cylinder_tangential_vectors(self):
        """
        Test the ``calc_cylinder_tangential_vectors`` method.
        """
        agrid = 1.
        offset = 0.5
        center = np.array(3 * [offset])
        node_indices = np.array([[0, 0, 0],
                                 [2, 0, 0],
                                 [0, 2, 0],
                                 [-2, 0, 0],
                                 [0, -2, 0]])
        ref_tangents = np.array([[0, 0, 0],
                                 [0, 1, 0],
                                 [-1, 0, 0],
                                 [0, -1, 0],
                                 [1, 0, 0]])
        tangents = espressomd.lb.calc_cylinder_tangential_vectors(
            center, agrid, offset, node_indices)
        np.testing.assert_array_almost_equal(tangents, ref_tangents)


if __name__ == "__main__":
    ut.main()
