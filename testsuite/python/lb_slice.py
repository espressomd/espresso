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


class LBTest:

    """This simple test first writes random numbers and then reads them
    to same slices of LB nodes and compares if the results are the same,
    shape-wise and value-wise.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.1
    np.random.seed(seed=42)

    def setUp(self):
        self.lb_fluid = self.lb_class(
            agrid=1., density=1., kinematic_viscosity=1.,
            tau=self.system.time_step, **self.lb_params)
        self.system.actors.add(self.lb_fluid)

    def tearDown(self):
        self.system.actors.clear()

    def test_slicing(self):
        lb_fluid = self.lb_fluid

        # array locked
        array = lb_fluid[1:-1:1, 5, 3:6].velocity
        with self.assertRaisesRegex(ValueError, "ESPResSo array properties return non-writable arrays"):
            array[0, 0, 0, 1] = 5.

        # density broadcast (with type conversion from int to double)
        lb_fluid[:, :, 0].density = 2
        np.testing.assert_array_almost_equal(
            np.copy(lb_fluid[:, :, 0].density), 2.)

        # velocity on test slice [:-1, :-1, -1]
        input_vel = np.random.rand(9, 9, 9, 3)
        lb_fluid[:-1, :-1, :-1].velocity = input_vel
        output_vel = lb_fluid[:-1, :-1, :-1].velocity
        np.testing.assert_array_almost_equal(input_vel, np.copy(output_vel))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of 'velocity' array \(9, 9, 9, 2\) does not match slice dimensions \(9, 9, 9, 3\)"):
            lb_fluid[:-1, :-1, :-1].velocity = input_vel[:, :, :, :2]

        # velocity broadcast (with type conversion from int to double)
        lb_fluid[:, :, 0].velocity = [1, 2, 3]
        np.testing.assert_array_almost_equal(
            np.copy(lb_fluid[:, :, 0].velocity), 10 * [10 * [[1, 2, 3]]])

        input_dens = np.random.rand(8, 3) + 1.
        lb_fluid[1:-1, 5, 3:6].density = input_dens
        output_dens = lb_fluid[1:-1, 5, 3:6].density
        np.testing.assert_array_almost_equal(np.copy(output_dens), input_dens)

        # population on test slice [:, :, :]
        input_pop = np.random.rand(10, 10, 10, 19)
        lb_fluid[:, :, :].population = input_pop
        output_pop = lb_fluid[:, :, :].population
        np.testing.assert_array_almost_equal(input_pop, np.copy(output_pop))

        with self.assertRaisesRegex(ValueError, r"Input-dimensions of 'population' array \(10, 10, 10, 5\) does not match slice dimensions \(10, 10, 10, 19\)"):
            lb_fluid[:, :, :].population = input_pop[:, :, :, :5]

        # pressure tensor on test slice [3, 6, 2:5]
        output_pressure_shape = lb_fluid[3, 6, 2:5].pressure_tensor.shape
        should_pressure_shape = (3, 3, 3)
        np.testing.assert_array_equal(
            output_pressure_shape, should_pressure_shape)

        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor' is read-only"):
            lb_fluid[3, 6, 2:5].pressure_tensor = np.zeros(
                should_pressure_shape)

        # pressure tensor non-equilibrium on test slice [3, 6, 2:5]
        output_pressure_shape = lb_fluid[3, 6, 2:5].pressure_tensor_neq.shape
        should_pressure_shape = (3, 3, 3)
        np.testing.assert_array_equal(
            output_pressure_shape, should_pressure_shape)

        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor_neq' is read-only"):
            lb_fluid[3, 6, 2:5].pressure_tensor_neq = np.zeros(
                should_pressure_shape)

        # boundary velocity on test slice [1:, 1:, 1:]
        output_boundary_shape = lb_fluid[1:, 1:, 1:].boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_equal(
            output_boundary_shape, should_boundary_shape)

        with self.assertRaisesRegex(TypeError, "Parameter 'values' must be an array_like of VelocityBounceBack or None"):
            lb_fluid[1:, 1:, 1:].boundary = np.zeros(should_boundary_shape)
        with self.assertRaisesRegex(TypeError, "Parameter 'values' must be an array_like of VelocityBounceBack or None"):
            lb_fluid[1:, 1:, 1:].boundary = np.array(
                [None, [1, 2, 3]], dtype=object)

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
            10 * [10 * [[1, 2, 3]]])

        # access out of bounds
        i = lb_fluid.shape[2] + 10
        lb_slice = lb_fluid[1, 2, i:i + 10]
        self.assertEqual(lb_slice.density.shape, (0,))
        self.assertIsInstance(lb_slice.density.dtype, object)
        with self.assertRaisesRegex(AttributeError, "Cannot set properties of an empty 'LBFluidSliceWalberla' object"):
            lb_slice.density = [1., 2., 3.]

    def test_iterator(self):
        lbslice_handle = self.lb_fluid[:, :, :]
        # arrange node indices using class methods
        lb_indices = [np.arange(self.lb_fluid.shape[i]) for i in range(3)]
        arranged_indices = list(itertools.product(*lb_indices))
        # arrange node indices using __iter__() enforced conversion
        iterator_indices = [node.index for node in lbslice_handle]
        # check the results correspond pairwise. order is implicitly preserved.
        np.testing.assert_array_equal(arranged_indices, iterator_indices)
        # use __eq()__ method form LBFluidRoutines()
        assert all([x == y for x, y in zip(
            arranged_indices, iterator_indices)])


@utx.skipIfMissingFeatures("WALBERLA")
class LBTestWalberlaDoublePrecisionCPU(LBTest, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {"single_precision": False}


@utx.skipIfMissingFeatures("WALBERLA")
class LBTestWalberlaSinglePrecisionCPU(LBTest, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {"single_precision": True}


if __name__ == "__main__":
    ut.main()
