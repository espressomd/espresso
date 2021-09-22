# Copyright (C) 2010-2019 The ESPResSo project
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
import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes
import unittest as ut
import unittest_decorators as utx
import numpy as np


@utx.skipIfMissingFeatures(["LB_BOUNDARIES", "LB_WALBERLA"])
class LBBoundaryVelocityTest(ut.TestCase):
    """
    Various tests to check the interaction of lb velocity boundary conditions and the fluid
    """

    lb_params = {'agrid': 0.6,
                 'dens': 0.5,
                 'visc': 3.2,
                 'tau': 0.7}
    system = espressomd.System(box_l=3 * [8 * lb_params['agrid']])
    system.time_step = lb_params['tau']
    system.cell_system.skin = 0.1

    lb_fluid = None

    def tearDown(self) -> None:
        self.system.actors.clear()
        self.system.lbboundaries.clear()

    def setUp(self) -> None:
        self.lb_fluid = espressomd.lb.LBFluidWalberla(**self.lb_params)
        self.system.actors.add(self.lb_fluid)

    def check_wall_slip(self, v_boundary):
        """
        Check that the fluid adopts the velocity set by the boundary conditions.
        """

        agrid = self.lb_params['agrid']
        wall_shape_left = espressomd.shapes.Wall(normal=[1, 2, 3], dist=0)
        wall_shape_right = espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-(self.system.box_l[0] - agrid))
        for shape in [wall_shape_left, wall_shape_right]:
            wall = espressomd.lbboundaries.LBBoundary(
                shape=shape, velocity=v_boundary)
            self.system.lbboundaries.add(wall)

        # velocity in front of the wall must be that of the boundary immediately
        self.system.integrator.run(100)

        for node in self.lb_fluid.nodes():
            if not node.is_boundary:
                print(node.velocity)

        v_fluid = self.lb_fluid[2, 0, 0].velocity
        np.testing.assert_allclose(v_fluid, v_boundary, atol=1E-5)

        # velocity in the middle needs to propagate first
        self.system.integrator.run(200)
        v_fluid = self.lb_fluid[2, 1, 3].velocity
        np.testing.assert_allclose(v_fluid, v_boundary, atol=1E-5)

    def test_wall_slip_parallel(self):
        v_boundary = [0, 0, 0.07]
        self.check_wall_slip(v_boundary)

    def test_wall_slip_nonparallel(self):
        v_boundary = [0.03, 0.05, 0.07]
        self.check_wall_slip(v_boundary)

    def test_boundary_readout(self):
        """
        Test the read part of the boundary property of lb nodes.
        """
        v_boundary = [0.03, 0.05, 0.07]
        wall_shape = espressomd.shapes.Wall(normal=[1, 0, 0], dist=self.lb_params['agrid'])
        wall = espressomd.lbboundaries.LBBoundary(
            shape=wall_shape, velocity=v_boundary)
        self.system.lbboundaries.add(wall)

        # check non_boundary node
        bound_cond = self.lb_fluid[4, 4, 4].boundary
        self.assertIsNone(bound_cond)

        # on boundary
        bound_cond = self.lb_fluid[0, 2, 4].boundary
        np.testing.assert_array_almost_equal(bound_cond.velocity, v_boundary)

        # TODO boundary_force

    def test_velocity_bounce_back_class(self):
        """
        Test setters and getters of :ref:`espressomd.lbboundaries.VelocityBounceBack`
        """
        with self.assertRaises(ValueError):
            bound_cond = espressomd.lbboundaries.VelocityBounceBack([1, 2, 3, 4])
        v = [1, 2, 17.4]
        bound_cond = espressomd.lbboundaries.VelocityBounceBack(v)
        np.testing.assert_array_almost_equal(bound_cond.velocity, v)

    def test_boundary_setting(self):
        """
        Test setting and un-setting individual lb boundary nodes.
        """
        v_boundary = [0.2, 0.1, 0.3]
        bound_cond = espressomd.lbboundaries.VelocityBounceBack(v_boundary)

        with self.assertRaises(ValueError):
            self.lb_fluid[1, 2, 3].boundary = 17

        self.lb_fluid[1, 2, 3].boundary = bound_cond
        np.testing.assert_array_almost_equal(self.lb_fluid[1, 2, 3].boundary.velocity, bound_cond.velocity)

        self.lb_fluid[1, 2, 3].boundary = None
        self.assertIsNone(self.lb_fluid[1, 2, 3].boundary)

    def test_nodes_in_shape(self):
        """
        Test if the get_nodes_in_shape method correctly identifies the grid points inside a cylinder
        """
        agrid = self.lb_params['agrid']
        cyl = espressomd.shapes.Cylinder(center=agrid * np.array([1.5, 2.5, 5.5]),
                                         axis=[0, 0, 1],
                                         length=2.1 * agrid,
                                         radius=0.5 * agrid)
        nodes_in_cyl = self.lb_fluid.get_nodes_in_shape(cyl)
        idxs_in_cyl = set(tuple(node.index) for node in nodes_in_cyl)

        idx_should_be = {(1, 2, 5), (1, 2, 4), (1, 2, 6)}
        self.assertSetEqual(idxs_in_cyl, idx_should_be)


if __name__ == "__main__":
    ut.main()
