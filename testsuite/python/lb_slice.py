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


@utx.skipIfMissingFeatures(["LB_BOUNDARIES"])
class LBSliceTest(ut.TestCase):

    """This simple test first writes random numbers and then reads them 
    to same slices of LB nodes and compares if the results are the same, 
    shape and value wise.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = .01
    system.cell_system.skin = 0.1

    def test_slicing(self):
        system = self.system

        lb_fluid = espressomd.lb.LBFluid(
            agrid=1.0, dens=1., visc=1., tau=0.01)
        system.actors.add(lb_fluid)

        "velocity on test slice [:-1, :-1, -1]"
        input_vel = np.random.rand(9, 9, 9, 3)
        lb_fluid[:-1, :-1, :-1].velocity = input_vel
        output_vel = lb_fluid[:-1, :-1, :-1].velocity
        np.testing.assert_array_almost_equal(input_vel, output_vel)

        "density on test slice [1:-1:2, 5, 3:6:2]"
        input_dens = np.random.rand(4, 1, 2)
        lb_fluid[1:-1:2, 5, 3:6:2].density = input_dens
        output_dens = lb_fluid[1:-1:2, 5, 3:6:2].density
        np.testing.assert_array_almost_equal(input_dens, output_dens)

        "population on test slice [:, :, :]"
        input_pop = np.random.rand(10, 10, 10, 19)
        lb_fluid[:, :, :].population = input_pop
        output_pop = lb_fluid[:, :, :].population
        np.testing.assert_array_almost_equal(input_pop, output_pop)

        "pressure tensor on test slice [3, 6, 2:5], should be of shape (1, 1, 3, 3, 3)"
        output_pressure_shape = lb_fluid[3, 6, 2:5].pressure_tensor.shape
        should_pressure_shape = (1, 1, 3, 3, 3)
        np.testing.assert_array_almost_equal(
            output_pressure_shape, should_pressure_shape)

        "pressure tensor neq on test slice [3, 6, 2:10], should be of shape (2, 1, 8, 3, 3)"
        output_pressure_neq_shape = lb_fluid[3:5,
                                             6:7, 2:10].pressure_tensor_neq.shape
        should_pressure_neq_shape = (2, 1, 8, 3, 3)
        np.testing.assert_array_almost_equal(
            output_pressure_neq_shape, should_pressure_neq_shape)

        "index on test slice [1, 1:5, 6:], should be of shape (1, 4, 4, 3)"
        output_index_shape = lb_fluid[1, 1:5, 6:].index.shape
        should_index_shape = (1, 4, 4, 3)
        np.testing.assert_array_almost_equal(
            output_index_shape, should_index_shape)

        "boundary on test slice [1:, 1:, 1:], should be of shape (9, 9, 9)"
        output_boundary_shape = lb_fluid[1:, 1:, 1:].boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_almost_equal(
            output_boundary_shape, should_boundary_shape)


if __name__ == "__main__":
    ut.main()
