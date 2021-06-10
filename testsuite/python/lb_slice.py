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
import numpy as np
import itertools


class LBSliceTest(ut.TestCase):

    """This simple test first writes random numbers and then reads them 
    to same slices of LB nodes and compares if the results are the same, 
    shape and value wise.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = .01
    system.cell_system.skin = 0.1
    np.random.seed(seed=42)
    lb_fluid = espressomd.lb.LBFluid(
        agrid=1.0, dens=1., visc=1., tau=0.01)
    system.actors.add(lb_fluid)

    def test_slicing(self):
        # array locked
        array = self.lb_fluid[1:-1:2, 5, 3:6:2].velocity
        with self.assertRaisesRegex(ValueError, "ESPResSo array properties return non-writable arrays"):
            array[0, 0, 0, 1] = 5.

        # velocity on test slice [:-1, :-1, -1]
        input_vel = np.random.rand(9, 9, 9, 3)
        self.lb_fluid[:-1, :-1, :-1].velocity = input_vel
        output_vel = self.lb_fluid[:-1, :-1, :-1].velocity
        np.testing.assert_array_almost_equal(input_vel, np.copy(output_vel))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of velocity array \(9, 9, 9, 2\) does not match slice dimensions \(9, 9, 9, 3\)"):
            self.lb_fluid[:-1, :-1, :-1].velocity = input_vel[:, :, :, :2]

        # velocity broadcast
        self.lb_fluid[:, :, 0].velocity = [1, 2, 3]
        np.testing.assert_array_almost_equal(
            np.copy(self.lb_fluid[:, :, 0].velocity), 10 * [10 * [[[1, 2, 3]]]])

        # density on test slice [1:-1:2, 5, 3:6:2]
        input_dens = np.random.rand(4, 1, 2)
        self.lb_fluid[1:-1:2, 5, 3:6:2].density = input_dens
        output_dens = self.lb_fluid[1:-1:2, 5, 3:6:2].density
        np.testing.assert_array_almost_equal(input_dens, np.copy(output_dens))

        # density broadcast
        self.lb_fluid[:, :, 0].density = 1.2
        np.testing.assert_array_almost_equal(
            np.copy(self.lb_fluid[:, :, 0].density), 1.2)

        # population on test slice [:, :, :]
        input_pop = np.random.rand(10, 10, 10, 19)
        self.lb_fluid[:, :, :].population = input_pop
        output_pop = self.lb_fluid[:, :, :].population
        np.testing.assert_array_almost_equal(input_pop, np.copy(output_pop))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of population array \(10, 10, 10, 5\) does not match slice dimensions \(10, 10, 10, 19\)"):
            self.lb_fluid[:, :, :].population = input_pop[:, :, :, :5]

        # pressure tensor on test slice [3, 6, 2:5]
        output_pressure_shape = self.lb_fluid[3, 6, 2:5].pressure_tensor.shape
        should_pressure_shape = (1, 1, 3, 3, 3)
        np.testing.assert_array_almost_equal(
            output_pressure_shape, should_pressure_shape)

        with self.assertRaises(NotImplementedError):
            self.lb_fluid[3, 6, 2:5].pressure_tensor = np.zeros(
                should_pressure_shape)

        # pressure tensor neq on test slice [3, 6, 2:10]
        output_pressure_neq_shape = self.lb_fluid[3:5,
                                                  6:7, 2:10].pressure_tensor_neq.shape
        should_pressure_neq_shape = (2, 1, 8, 3, 3)
        np.testing.assert_array_almost_equal(
            output_pressure_neq_shape, should_pressure_neq_shape)

        with self.assertRaises(NotImplementedError):
            self.lb_fluid[3:5, 6:7, 2:10].pressure_tensor_neq = np.zeros(
                output_pressure_neq_shape)

        # index on test slice [1, 1:5, 6:]
        output_index_shape = self.lb_fluid[1, 1:5, 6:].index.shape
        should_index_shape = (1, 4, 4, 3)
        np.testing.assert_array_almost_equal(
            output_index_shape, should_index_shape)

        with self.assertRaisesRegex(AttributeError, "attribute 'index' of 'espressomd.lb.LBFluidRoutines' objects is not writable"):
            self.lb_fluid[1, 1:5, 6:].index = np.zeros(output_index_shape)

        # boundary on test slice [1:, 1:, 1:]
        if espressomd.has_features('LB_BOUNDARIES'):
            output_boundary_shape = self.lb_fluid[1:, 1:, 1:].boundary.shape
            should_boundary_shape = (9, 9, 9)
            np.testing.assert_array_almost_equal(
                output_boundary_shape, should_boundary_shape)

            with self.assertRaises(NotImplementedError):
                self.lb_fluid[1:, 1:, 1:].boundary = np.zeros(
                    should_boundary_shape)

    def test_iterator(self):
        lbslice_handle = self.lb_fluid[:, :, :]
        # arrange node indices using class methods
        i_handle, j_handle, k_handle = lbslice_handle.x_indices, lbslice_handle.y_indices, lbslice_handle.z_indices
        arranged_indices = [
            (x, y, z) for (
                x, y, z) in itertools.product(
                i_handle, j_handle, k_handle)]
        # arrange node indices using __iter__() enforced converstion
        iterator_indices = [x.index for x in lbslice_handle]
        # check the results correspond pairwise. order implicitly preserved.
        # uses __eq()__ method form LBFluidRoutines()
        assert all([x == y for x, y in zip(
            arranged_indices, iterator_indices)])


if __name__ == "__main__":
    ut.main()
