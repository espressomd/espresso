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
    system.part.add(pos=[0, 0, 0])

    def setUp(self):
        self.system.box_l = [12.0, 12.0, 12.0]

    def tearDown(self):
        self.system.actors.clear()

    def assert_copy_is_writable(self, array):
        cpy = np.copy(array)
        self.assertTrue(cpy.flags.writeable)

    def test_common(self):
        self.assert_operator_usage_raises(self.system.part[0].pos)
        self.assert_operator_usage_raises(self.system.part[0].v)
        self.assert_operator_usage_raises(self.system.part[0].f)
        self.assert_operator_usage_raises(self.system.part[0].pos_folded)

        self.assert_operator_usage_raises(self.system.box_l)

        check_array_writable(self.system.part[0].pos)
        check_array_writable(self.system.part[0].v)
        check_array_writable(self.system.part[0].f)
        check_array_writable(self.system.box_l)

        self.assert_copy_is_writable(self.system.part[0].pos)
        self.assert_copy_is_writable(self.system.part[0].v)
        self.assert_copy_is_writable(self.system.part[0].f)
        self.assert_copy_is_writable(self.system.part[0].pos_folded)

        self.assert_copy_is_writable(self.system.box_l)

    @utx.skipIfMissingFeatures(["ROTATION"])
    def test_rotation(self):
        self.assert_operator_usage_raises(self.system.part[0].omega_lab)
        self.assert_operator_usage_raises(self.system.part[0].quat)
        self.assert_operator_usage_raises(self.system.part[0].rotation)
        self.assert_operator_usage_raises(self.system.part[0].omega_body)
        self.assert_operator_usage_raises(self.system.part[0].torque_lab)
        if espressomd.has_features("EXTERNAL_FORCES"):
            self.assert_operator_usage_raises(self.system.part[0].ext_torque)

        check_array_writable(self.system.part[0].quat)
        check_array_writable(self.system.part[0].omega_lab)
        check_array_writable(self.system.part[0].rotation)
        check_array_writable(self.system.part[0].omega_body)
        check_array_writable(self.system.part[0].torque_lab)

        if espressomd.has_features("EXTERNAL_FORCES"):
            check_array_writable(self.system.part[0].ext_torque)

        self.assert_copy_is_writable(self.system.part[0].omega_lab)
        self.assert_copy_is_writable(self.system.part[0].quat)
        self.assert_copy_is_writable(self.system.part[0].rotation)
        self.assert_copy_is_writable(self.system.part[0].omega_body)
        self.assert_copy_is_writable(self.system.part[0].torque_lab)
        if espressomd.has_features("EXTERNAL_FORCES"):
            self.assert_copy_is_writable(self.system.part[0].ext_torque)

    @utx.skipIfMissingFeatures(["ROTATIONAL_INERTIA"])
    def test_rotational_inertia(self):
        self.assert_operator_usage_raises(self.system.part[0].rinertia)
        check_array_writable(self.system.part[0].rinertia)
        self.assert_copy_is_writable(self.system.part[0].rinertia)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_external_forces(self):
        self.assert_operator_usage_raises(self.system.part[0].ext_force)
        self.assert_operator_usage_raises(self.system.part[0].fix)

        check_array_writable(self.system.part[0].ext_force)
        check_array_writable(self.system.part[0].fix)

        self.assert_copy_is_writable(self.system.part[0].ext_force)
        self.assert_copy_is_writable(self.system.part[0].fix)

    @utx.skipIfMissingFeatures(["ROTATION", "PARTICLE_ANISOTROPY"])
    def test_rot_aniso(self):
        self.assert_operator_usage_raises(self.system.part[0].gamma_rot)

        check_array_writable(self.system.part[0].gamma_rot)

        self.assert_copy_is_writable(self.system.part[0].gamma_rot)

    def test_lb(self):
        lbf = espressomd.lb.LBFluid(agrid=0.5, dens=1, visc=1, tau=0.01)
        self.system.actors.add(lbf)

        self.assert_operator_usage_raises(lbf[0, 0, 0].velocity)
        self.assert_operator_usage_raises(lbf[0, 0, 0].stress)
        self.assert_operator_usage_raises(lbf[0, 0, 0].stress_neq)
        self.assert_operator_usage_raises(lbf[0, 0, 0].population)

    @utx.skipIfMissingFeatures(["LANGEVIN_PER_PARTICLE",
                                "PARTICLE_ANISOTROPY"])
    def test_langevinpp_aniso(self):
        self.assert_operator_usage_raises(self.system.part[0].gamma)

        check_array_writable(self.system.part[0].gamma)

        self.assert_copy_is_writable(self.system.part[0].gamma)

    @utx.skipIfMissingFeatures(["DIPOLES"])
    def test_dipoles(self):
        self.assert_operator_usage_raises(self.system.part[0].dip)

        check_array_writable(self.system.part[0].dip)

        self.assert_copy_is_writable(self.system.part[0].dip)

    @utx.skipIfMissingFeatures(["EXCLUSIONS"])
    def test_exclusions(self):
        self.assert_operator_usage_raises(self.system.part[0].exclusions)

    def test_partial_periodic(self):
        self.assert_operator_usage_raises(self.system.periodicity)

        check_array_writable(self.system.periodicity)

        self.assert_copy_is_writable(self.system.periodicity)


if __name__ == "__main__":
    ut.main()
