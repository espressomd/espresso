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
import unittest as ut
import unittest_decorators as utx
import numpy as np
import itertools


@utx.skipIfMissingFeatures("WALBERLA")
class LBSliceTest(ut.TestCase):

    """This simple test first writes random numbers and then reads them
    to same slices of LB nodes and compares if the results are the same,
    shape-wise and value-wise.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = .01
    system.cell_system.skin = 0.1
    np.random.seed(seed=42)

    @classmethod
    def setUpClass(cls):
        cls.lb_fluid = espressomd.lb.LBFluidWalberla(
            agrid=1.0, density=1., viscosity=1., tau=0.01)
        cls.system.actors.add(cls.lb_fluid)

    def test_slicing(self):
        lb_fluid = self.lb_fluid

        # array locked
        array = lb_fluid[1:-1:2, 5, 3:6:2].velocity
        with self.assertRaisesRegex(ValueError, "ESPResSo array properties return non-writable arrays"):
            array[0, 0, 0, 1] = 5.

        # velocity on test slice [:-1, :-1, -1]
        input_vel = np.random.rand(9, 9, 9, 3)
        lb_fluid[:-1, :-1, :-1].velocity = input_vel
        output_vel = lb_fluid[:-1, :-1, :-1].velocity
        np.testing.assert_array_almost_equal(input_vel, np.copy(output_vel))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of 'velocity' array \(9, 9, 9, 2\) does not match slice dimensions \(9, 9, 9, 3\)"):
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

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of 'population' array \(10, 10, 10, 5\) does not match slice dimensions \(10, 10, 10, 19\)"):
            lb_fluid[:, :, :].population = input_pop[:, :, :, :5]

        # pressure tensor on test slice [3, 6, 2:5]
        output_pressure_shape = lb_fluid[3, 6, 2:5].pressure_tensor.shape
        should_pressure_shape = (1, 1, 3, 3, 3)
        np.testing.assert_array_equal(
            output_pressure_shape, should_pressure_shape)

        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor' is read-only"):
            lb_fluid[3, 6, 2:5].pressure_tensor = np.zeros(
                should_pressure_shape)

        # index on test slice [1, 1:5, 6:]
        output_index_shape = lb_fluid[1, 1:5, 6:].index.shape
        should_index_shape = (1, 4, 4, 3)
        np.testing.assert_array_equal(
            output_index_shape, should_index_shape)

        with self.assertRaisesRegex(RuntimeError, "Parameter 'index' is read-only"):
            lb_fluid[1, 1:5, 6:].index = np.zeros(output_index_shape)

        # boundary velocity on test slice [1:, 1:, 1:]
        output_boundary_shape = lb_fluid[1:, 1:, 1:].boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_equal(
            output_boundary_shape, should_boundary_shape)

        with self.assertRaisesRegex(TypeError, "values must be instances of VelocityBounceBack or None"):
            lb_fluid[1:, 1:, 1:].boundary = np.zeros(should_boundary_shape)

        vbb_ref = espressomd.lb.VelocityBounceBack([1e-6, 2e-6, 3e-6])
        lb_fluid[1:2, 1:, 0].boundary = vbb_ref
        lb_fluid[1:2, 2:, 0].boundary = None
        for vbb in lb_fluid[1:2, 1, 0].boundary.flatten():
            np.testing.assert_array_almost_equal(
                vbb.velocity, vbb_ref.velocity)
        for vbb in lb_fluid[1:2, 2, 0:2].boundary.flatten():
            self.assertIsNone(vbb)

        # is_boundary on test slice [1:, 1:, 1:]
        output_boundary_shape = lb_fluid[1:, 1:, 1:].is_boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_equal(
            output_boundary_shape, should_boundary_shape)

        with self.assertRaisesRegex(RuntimeError, "Property 'is_boundary' is read-only"):
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

        # access out of bounds
        i = lb_fluid.shape[2] + 10
        lb_slice = lb_fluid[1, 2, i:i + 10]
        self.assertEqual(lb_slice.density.shape, (0,))
        self.assertIsInstance(lb_slice.density.dtype, object)
        self.assertEqual(list(lb_slice.indices[0]), [1])
        self.assertEqual(list(lb_slice.indices[1]), [2])
        self.assertEqual(list(lb_slice.indices[2]), [])
        np.testing.assert_array_equal(
            np.copy(lb_fluid[1, 2, 8:i].index),
            np.copy(lb_fluid[1, 2, 8:].index))
        with self.assertRaisesRegex(AttributeError, "Cannot set properties of an empty 'LBFluidSliceWalberla' object"):
            lb_slice.density = [1., 2., 3.]

    def test_iterator(self):
        lbslice_handle = self.lb_fluid[:, :, :]
        # arrange node indices using class methods
        arranged_indices = list(itertools.product(*lbslice_handle.indices))
        # arrange node indices using __iter__() enforced conversion
        iterator_indices = [node.index for node in lbslice_handle]
        # check the results correspond pairwise. order is implicitly preserved.
        np.testing.assert_array_equal(arranged_indices, iterator_indices)
        # use __eq()__ method form LBFluidRoutines()
        assert all([x == y for x, y in zip(
            arranged_indices, iterator_indices)])


if __name__ == "__main__":
    ut.main()
