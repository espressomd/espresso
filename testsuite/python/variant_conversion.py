# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import sys
import unittest as ut
import numpy as np

from espressomd import script_interface


@script_interface.script_interface_register
class VariantTester(script_interface.ScriptInterfaceHelper):
    _so_name = "Testing::VariantTester"
    _so_creation_policy = "LOCAL"


class test_variant_conversion(ut.TestCase):

    """
    Test the variant conversion functions.
    This needs the c++ helper class ScriptInterface::Testing::VariantTester.
    """
    vt = VariantTester()

    def check_type_and_value(self, stype, svalue, value):
        self.assertEqual(type(value), stype)
        self.assertEqual(value, svalue)

    def test_bool_return(self):
        """Check that bool return values work corrcetly,
        because they are needed for the other tests."""
        self.assertTrue(self.vt.call_method("true"))
        self.assertFalse(self.vt.call_method("false"))

    def test_mirror(self):
        for v in 1.1, 0.6, 0.1232, 1, [1, 1.1], "blabla":
            self.assertEqual(v, self.vt.call_method("mirror", value=v))

    def test_default(self):
        """ Check that a default constructed Variant translates to None.
        """
        self.assertTrue(self.vt.call_method("default") is None)

    def test_flat(self):
        ret = self.vt.call_method("flat")
        self.assertTrue(isinstance(ret, list))
        self.check_type_and_value(bool, True, ret[0])
        self.check_type_and_value(str, 'a string', ret[1])
        self.check_type_and_value(float, 3.14159, ret[2])
        self.check_type_and_value(list, [3, 1, 4, 1, 5], ret[3])
        self.check_type_and_value(list, [1.1, 2.2, 3.3], ret[4])
        # Empty object (ObjectId()) should map to None
        self.assertEqual(ret[5], None)
        # An actual object
        self.assertTrue(isinstance(ret[6], script_interface.PScriptInterface))
        self.assertEqual(ret[6].name(), "Testing::VariantTester")

    def test_recu(self):
        ret = self.vt.call_method("recursive", max_level=5)
        self.check_type_and_value(str, 'end', ret[1][1][1][1][1])

    def test_mixed(self):
        ret = self.vt.call_method("mixed")
        self.check_type_and_value(int, 1, ret[0])
        self.check_type_and_value(str, 'another string', ret[1])
        self.check_type_and_value(str, 'end', ret[2][1][1][1][1][1])

    def test_parameter_types(self):
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="none", value=None))
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="bool", value=True))
        self.assertTrue(self.vt.call_method(
            "check_parameter_type", type="int", value=42))
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="int", value=np.int64(42)))
        self.assertTrue(self.vt.call_method(
            "check_parameter_type", type="string", value='blub'))
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="double", value=12.5))
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="double", value=np.float64(12.5)))
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="objectid", value=self.vt))
        self.assertTrue(self.vt.call_method(
            "check_parameter_type", type="double_vector", value=[1.1, 2.2, 3.3]))
        self.assertTrue(
            self.vt.call_method("check_parameter_type", type="int_vector", value=[1, 2, 3]))
        self.assertTrue(self.vt.call_method(
            "check_parameter_type", type="vector", value=[1, 'string', True]))

if __name__ == "__main__":
    ut.main()
