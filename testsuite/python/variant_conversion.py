from __future__ import print_function
import sys
import unittest as ut

from espressomd import script_interface

class test_variant_conversion(ut.TestCase):
    """
    Test the variant conversion functions.
    This needs the c++ helper class ScriptInterface::Testing::VariantTester.
    """
    vt = script_interface.PScriptInterface("Testing::VariantTester")

    def check_type_and_value(self, stype, svalue, value):
        assert(type(value) == stype)
        assert(value == svalue)

    def test_bool_return(self):
        """Check that bool return values work corrcetly,
        because they are needed for the other tests."""
        assert(self.vt.call_method("true") == True)
        assert(self.vt.call_method("false") == False)

    def test_flat(self):
        ret = self.vt.call_method("flat")
        assert(isinstance(ret, list))
        self.check_type_and_value(bool, True, ret[0])
        self.check_type_and_value(str, 'a string', ret[1])
        self.check_type_and_value(float, 3.14159, ret[2])
        self.check_type_and_value(list, [3,1,4,1,5], ret[3])
        self.check_type_and_value(list, [1.1,2.2,3.3], ret[4])
        # Empty object (ObjectId()) should map to None
        assert(ret[5] == None)
        # An actual object
        assert(isinstance(ret[6], script_interface.PScriptInterface))
        assert(ret[6].name() == "Testing::VariantTester")

    def test_recu(self):
        ret = self.vt.call_method("recursive", max_level=5)
        self.check_type_and_value(str, 'end', ret[1][1][1][1][1])

    def test_mixed(self):
        ret = self.vt.call_method("mixed")
        self.check_type_and_value(int, 1, ret[0])
        self.check_type_and_value(str, 'another string', ret[1])
        self.check_type_and_value(str, 'end', ret[2][1][1][1][1][1])

    def test_parameter_types(self):
        assert(self.vt.call_method("check_parameter_type", type="bool", value=True))
        assert(self.vt.call_method("check_parameter_type", type="int", value=42))
        assert(self.vt.call_method("check_parameter_type", type="string", value='blub'))
        assert(self.vt.call_method("check_parameter_type", type="double", value=12.5))
        assert(self.vt.call_method("check_parameter_type", type="objectid", value=self.vt))
        assert(self.vt.call_method("check_parameter_type", type="double_vector", value=[1.1, 2.2, 3.3]))
        assert(self.vt.call_method("check_parameter_type", type="int_vector", value=[1,2,3]))
        assert(self.vt.call_method("check_parameter_type", type="vector", value=[1,'string',True]))

if __name__ == "__main__":
    ut.main()
