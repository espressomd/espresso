#
# Copyright (C) 2020 The ESPResSo project
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

import espressomd.utils as utils


class espresso_utils(ut.TestCase):

    def test_check_type_or_throw_except(self):
        with self.assertRaisesRegex(
                ValueError, 'A -- 4 values were given but 3 were expected.'):
            utils.check_type_or_throw_except([1, 2, 3, 4], 3, float, 'A')
        with self.assertRaisesRegex(
                ValueError, 'B -- 2 values were given but 3 were expected.'):
            utils.check_type_or_throw_except([1, 2], 3, float, 'B')
        with self.assertRaisesRegex(
                ValueError, 'C -- A single value was given but 3 were expected.'):
            utils.check_type_or_throw_except(1, 3, float, 'C')
        with self.assertRaisesRegex(
                ValueError, 'D -- Item 1 was of type str'):
            utils.check_type_or_throw_except([1, '2', '3'], 3, float, 'D')
        # the following statements should not raise any exception
        try:
            utils.check_type_or_throw_except([1, 2, 3], 3, float, '')
            utils.check_type_or_throw_except(np.array([1, 2]), 2, float, '')
            utils.check_type_or_throw_except(np.array(2 * [True]), 2, bool, '')
            utils.check_type_or_throw_except(np.array([1, 2])[0], 1, float, '')
            utils.check_type_or_throw_except(np.array([True])[0], 1, bool, '')
            utils.check_type_or_throw_except(np.array(['12'])[0], 1, str, '')
        except ValueError as err:
            self.fail(f'check_type_or_throw_except raised ValueError("{err}")')

    def test_is_valid_type(self):
        # basic types
        self.assertFalse(utils.is_valid_type(None, int))
        self.assertFalse(utils.is_valid_type('12', int))
        self.assertFalse(utils.is_valid_type(0.99, int))
        self.assertFalse(utils.is_valid_type(12, float))
        self.assertFalse(utils.is_valid_type(1234, str))
        self.assertTrue(utils.is_valid_type(1.0, float))
        self.assertTrue(utils.is_valid_type(12345, int))
        self.assertTrue(utils.is_valid_type('123', str))
        self.assertTrue(utils.is_valid_type(np.array([123.])[0], float))
        self.assertTrue(utils.is_valid_type(np.array([1234])[0], int))
        self.assertTrue(utils.is_valid_type(np.array([True])[0], bool))
        # numpy types
        self.assertTrue(utils.is_valid_type(
            np.array([12], dtype=int)[0], int))
        self.assertTrue(utils.is_valid_type(
            np.array([12], dtype=int)[0], int))
        self.assertTrue(utils.is_valid_type(
            np.array([1.], dtype=float)[0], float))
        self.assertTrue(utils.is_valid_type(
            np.array([1.], dtype=np.float64)[0], float))

    def test_nesting_level(self):
        self.assertEqual(utils.nesting_level(12345), 0)
        self.assertEqual(utils.nesting_level('123'), 0)
        self.assertEqual(utils.nesting_level((1, )), 1)
        self.assertEqual(utils.nesting_level([1, ]), 1)
        self.assertEqual(utils.nesting_level([[1]]), 2)


if __name__ == "__main__":
    ut.main()
