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

import unittest as ut

import numpy as np

import espressomd
import espressomd.lb
import espressomd.utils
import unittest_decorators as utx


class ArrayCommon(ut.TestCase):
    def assert_operator_usage_raises(self, array):
        with self.assertRaises(ValueError):
            array[0] = 0
        with self.assertRaises(ValueError):
            array += [1, 1, 1]
        with self.assertRaises(ValueError):
            array -= [1, 1, 1]
        with self.assertRaises(ValueError):
            array *= [1, 1, 1]
        with self.assertRaises(ValueError):
            array /= [1, 1, 1]
        with self.assertRaises(ValueError):
            array //= [1, 1, 1]
        with self.assertRaises(ValueError):
            array %= [1, 1, 1]
        with self.assertRaises(ValueError):
            array **= [1, 1, 1]
        with self.assertRaises(ValueError):
            array <<= [1, 1, 1]
        with self.assertRaises(ValueError):
            array >>= [1, 1, 1]
        with self.assertRaises(ValueError):
            array &= [1, 1, 1]
        with self.assertRaises(ValueError):
            array |= [1, 1, 1]
        with self.assertRaises(ValueError):
            array ^= [1, 1, 1]


class ArrayLockedTest(ArrayCommon):
    def test_locked_operators(self):
        array = espressomd.utils.array_locked([1., 2., 3.])
        self.assert_operator_usage_raises(array)

    def test_unlocked_operators(self):
        array = espressomd.utils.array_locked([1, 2, 3])
        array2 = espressomd.utils.array_locked([4, 5, 6])
        add = array + array2
        sub = array - array2

        self.assertIsInstance(add, np.ndarray)
        self.assertIsInstance(sub, np.ndarray)

        self.assertTrue(add.flags.writeable)
        self.assertTrue(sub.flags.writeable)

        np.testing.assert_array_equal(
            add, np.add(np.copy(array), np.copy(array2)))
        np.testing.assert_array_equal(
            sub, np.subtract(
                np.copy(array), np.copy(array2)))
        np.testing.assert_array_equal(sub, -(array2 - array))

    def test_copy_is_writeable(self):
        array = np.copy(espressomd.utils.array_locked([1, 2, 3]))
        self.assertTrue(array.flags.writeable)

    def test_setter(self):
        array = espressomd.utils.array_locked([1, 2, 3])
        array = [4, 5, 6]
        np.testing.assert_array_equal(array, [4, 5, 6])


def check_array_writable(array):
    value = np.random.random(array.shape[0]).astype(type(array[0]))
    array = value
    np.testing.assert_array_almost_equal(np.copy(array), value)


class ArrayPropertyTest(ArrayCommon):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.box_l = [12.0, 12.0, 12.0]
    system.time_step = 0.01
    system.cell_system.skin = 0.01
    partcl = system.part.add(pos=[0, 0, 0])

    def setUp(self):
        self.system.box_l = [12.0, 12.0, 12.0]

    def tearDown(self):
        self.system.actors.clear()

    def assert_copy_is_writable(self, array):
        cpy = np.copy(array)
        self.assertTrue(cpy.flags.writeable)

    def test_common(self):
        self.assert_operator_usage_raises(self.partcl.pos)
        self.assert_operator_usage_raises(self.partcl.v)
        self.assert_operator_usage_raises(self.partcl.f)
        self.assert_operator_usage_raises(self.partcl.pos_folded)

        self.assert_operator_usage_raises(self.system.box_l)

        check_array_writable(self.partcl.pos)
        check_array_writable(self.partcl.v)
        check_array_writable(self.partcl.f)
        check_array_writable(self.system.box_l)

        self.assert_copy_is_writable(self.partcl.pos)
        self.assert_copy_is_writable(self.partcl.v)
        self.assert_copy_is_writable(self.partcl.f)
        self.assert_copy_is_writable(self.partcl.pos_folded)

        self.assert_copy_is_writable(self.system.box_l)

    @utx.skipIfMissingFeatures(["ROTATION"])
    def test_rotation(self):
        self.assert_operator_usage_raises(self.partcl.omega_lab)
        self.assert_operator_usage_raises(self.partcl.quat)
        self.assert_operator_usage_raises(self.partcl.rotation)
        self.assert_operator_usage_raises(self.partcl.omega_body)
        self.assert_operator_usage_raises(self.partcl.torque_lab)
        if espressomd.has_features("EXTERNAL_FORCES"):
            self.assert_operator_usage_raises(self.partcl.ext_torque)

        check_array_writable(self.partcl.quat)
        check_array_writable(self.partcl.omega_lab)
        check_array_writable(self.partcl.rotation)
        check_array_writable(self.partcl.omega_body)
        check_array_writable(self.partcl.torque_lab)

        if espressomd.has_features("EXTERNAL_FORCES"):
            check_array_writable(self.partcl.ext_torque)

        self.assert_copy_is_writable(self.partcl.omega_lab)
        self.assert_copy_is_writable(self.partcl.quat)
        self.assert_copy_is_writable(self.partcl.rotation)
        self.assert_copy_is_writable(self.partcl.omega_body)
        self.assert_copy_is_writable(self.partcl.torque_lab)
        if espressomd.has_features("EXTERNAL_FORCES"):
            self.assert_copy_is_writable(self.partcl.ext_torque)

    @utx.skipIfMissingFeatures(["ROTATIONAL_INERTIA"])
    def test_rotational_inertia(self):
        self.assert_operator_usage_raises(self.partcl.rinertia)
        check_array_writable(self.partcl.rinertia)
        self.assert_copy_is_writable(self.partcl.rinertia)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_external_forces(self):
        self.assert_operator_usage_raises(self.partcl.ext_force)
        self.assert_operator_usage_raises(self.partcl.fix)

        check_array_writable(self.partcl.ext_force)
        check_array_writable(self.partcl.fix)

        self.assert_copy_is_writable(self.partcl.ext_force)
        self.assert_copy_is_writable(self.partcl.fix)

    @utx.skipIfMissingFeatures(["ROTATION", "THERMOSTAT_PER_PARTICLE",
                                "PARTICLE_ANISOTROPY"])
    def test_rot_aniso(self):
        self.assert_operator_usage_raises(self.partcl.gamma_rot)

        check_array_writable(self.partcl.gamma_rot)

        self.assert_copy_is_writable(self.partcl.gamma_rot)

    @utx.skipIfMissingFeatures("WALBERLA")
    def test_lb(self):
        lbf = espressomd.lb.LBFluidWalberla(
            agrid=0.5, density=1., kinematic_viscosity=1., tau=0.01)
        self.system.actors.add(lbf)

        self.assert_operator_usage_raises(lbf[0, 0, 0].velocity)
        self.assert_operator_usage_raises(lbf[0, 0, 0].pressure_tensor)
        self.assert_operator_usage_raises(lbf[0, 0, 0].population)

    @utx.skipIfMissingFeatures(["THERMOSTAT_PER_PARTICLE",
                                "PARTICLE_ANISOTROPY"])
    def test_thermostat_per_particle_aniso(self):
        self.assert_operator_usage_raises(self.partcl.gamma)

        check_array_writable(self.partcl.gamma)

        self.assert_copy_is_writable(self.partcl.gamma)

    @utx.skipIfMissingFeatures(["DIPOLES"])
    def test_dipoles(self):
        self.assert_operator_usage_raises(self.partcl.dip)

        check_array_writable(self.partcl.dip)

        self.assert_copy_is_writable(self.partcl.dip)

    @utx.skipIfMissingFeatures(["EXCLUSIONS"])
    def test_exclusions(self):
        self.assert_operator_usage_raises(self.partcl.exclusions)

    def test_partial_periodic(self):
        self.assert_operator_usage_raises(self.system.periodicity)

        check_array_writable(self.system.periodicity)

        self.assert_copy_is_writable(self.system.periodicity)


if __name__ == "__main__":
    ut.main()
