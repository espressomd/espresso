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
import unittest as ut
import unittest_decorators as utx
import numpy as np


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBSliceTest(ut.TestCase):

    """This simple test first writes random numbers and then reads them 
    to same slices of LB nodes and compares if the results are the same, 
    shape and value wise.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = .01
    system.cell_system.skin = 0.1
    np.random.seed(seed=42)

    def test_slicing(self):
        system = self.system

        lb_fluid = espressomd.lb.LBFluidWalberla(
            agrid=1.0, dens=1., visc=1., tau=0.01)
        system.actors.add(lb_fluid)

        # array locked
        array = lb_fluid[1:-1:2, 5, 3:6:2].velocity
        with self.assertRaisesRegex(ValueError, "ESPResSo array properties return non-writable arrays"):
            array[0, 0, 0, 1] = 5.

        # velocity on test slice [:-1, :-1, -1]
        input_vel = np.random.rand(9, 9, 9, 3)
        lb_fluid[:-1, :-1, :-1].velocity = input_vel
        output_vel = lb_fluid[:-1, :-1, :-1].velocity
        np.testing.assert_array_almost_equal(input_vel, np.copy(output_vel))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of velocity array \(9, 9, 9, 2\) does not match slice dimensions \(9, 9, 9, 3\)"):
            lb_fluid[:-1, :-1, :-1].velocity = input_vel[:, :, :, :2]

        # velocity broadcast
        lb_fluid[:, :, 0].velocity = [1, 2, 3]
        np.testing.assert_array_almost_equal(
            np.copy(lb_fluid[:, :, 0].velocity), 10 * [10 * [[[1, 2, 3]]]])

        # density on test slice [1:-1:2, 5, 3:6:2]
        input_dens = np.random.rand(4, 1, 2)
        lb_fluid[1:-1:2, 5, 3:6:2].density = input_dens
        output_dens = lb_fluid[1:-1:2, 5, 3:6:2].density
        np.testing.assert_array_almost_equal(input_dens, np.copy(output_dens))

        # density broadcast
        lb_fluid[:, :, 0].density = 1.2
        np.testing.assert_array_almost_equal(
            np.copy(lb_fluid[:, :, 0].density), 1.2)

        # population on test slice [:, :, :]
        input_pop = np.random.rand(10, 10, 10, 19)
        lb_fluid[:, :, :].population = input_pop
        output_pop = lb_fluid[:, :, :].population
        np.testing.assert_array_almost_equal(input_pop, np.copy(output_pop))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of population array \(10, 10, 10, 5\) does not match slice dimensions \(10, 10, 10, 19\)"):
            lb_fluid[:, :, :].population = input_pop[:, :, :, :5]

        # TODO Walberla: uncomment this block when pressure tensor is available
#        # pressure tensor on test slice [3, 6, 2:5]
#        output_pressure_shape = lb_fluid[3, 6, 2:5].pressure_tensor.shape
#        should_pressure_shape = (1, 1, 3, 3, 3)
#        np.testing.assert_array_almost_equal(
#            output_pressure_shape, should_pressure_shape)

#        with self.assertRaises(NotImplementedError):
#            lb_fluid[3, 6, 2:5].pressure_tensor = np.zeros(
#                should_pressure_shape)

        # index on test slice [1, 1:5, 6:]
        output_index_shape = lb_fluid[1, 1:5, 6:].index.shape
        should_index_shape = (1, 4, 4, 3)
        np.testing.assert_array_almost_equal(
            output_index_shape, should_index_shape)

        with self.assertRaisesRegex(AttributeError, "attribute 'index' of 'espressomd.lb.LBFluidRoutines' objects is not writable"):
            lb_fluid[1, 1:5, 6:].index = np.zeros(output_index_shape)

        # is_boundary on test slice [1:, 1:, 1:]
        output_boundary_shape = lb_fluid[1:, 1:, 1:].is_boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_almost_equal(
            output_boundary_shape, should_boundary_shape)

        with self.assertRaises(NotImplementedError):
            lb_fluid[1:, 1:, 1:].is_boundary = np.zeros(should_boundary_shape)

        # last_applied_force on test slice [:-1, :-1, -1]
        input_laf = np.random.rand(9, 9, 9, 3)
        lb_fluid[:-1, :-1, :-1].last_applied_force = input_laf
        output_laf = lb_fluid[:-1, :-1, :-1].last_applied_force
        np.testing.assert_array_almost_equal(input_laf, np.copy(output_laf))

        # last_applied_force broadcast
        lb_fluid[:, :, 0].last_applied_force = [1, 2, 3]
        np.testing.assert_array_almost_equal(
            np.copy(lb_fluid[:, :, 0].last_applied_force),
            10 * [10 * [[[1, 2, 3]]]])


if __name__ == "__main__":
    ut.main()
