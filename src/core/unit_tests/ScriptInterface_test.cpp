/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <type_traits>

#define BOOST_TEST_MODULE ScriptInterface test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/ScriptInterface.hpp"

using std::string;
using std::map;
using std::vector;

using namespace ScriptInterface;

namespace Testing {

/**
 * @brief Mock to test ScriptInterface.
 */
struct ScriptInterfaceTest : public ScriptInterface::ScriptInterfaceBase {
  map<string, Variant> get_parameters() const override {
    map<string, Variant> ret;

    ret["bool_opt"] = bool_opt;
    ret["integer"] = integer;
    ret["bool_req"] = bool_req;
    ret["vec_double"] = vec_double;
    ret["vec_int"] = vec_int;

    return ret;
  }

  /* Not needed for testing */
  map<string, Parameter> valid_parameters() const override { return {}; }

  void set_parameter(const string &name, const Variant &value) override {
    SET_PARAMETER_HELPER("bool_opt", bool_opt);
    SET_PARAMETER_HELPER("integer", integer);
    SET_PARAMETER_HELPER("bool_req", bool_req);
    SET_PARAMETER_HELPER("vec_double", vec_double);
    SET_PARAMETER_HELPER("vec_int", vec_int);
  }

  Variant call_method(const std::string &name,
                      const VariantMap &params) override {
    if (name == "test_method") {
      method_called = true;
    }

    return true;
  }

  bool method_called{false};
  bool bool_opt, bool_req;
  int integer;
  vector<double> vec_double;
  vector<int> vec_int;
};
} /* namespace Testing */

using namespace Testing;

BOOST_AUTO_TEST_CASE(non_copyable) {
  static_assert(!std::is_copy_constructible<ScriptInterfaceBase>::value, "");
  static_assert(!std::is_copy_assignable<ScriptInterfaceBase>::value, "");
}

/**
 * We check the default implementations of set_parameters
 * and get_parameter of ScriptInterface (this is the only
 * logic in the class).
 */
BOOST_AUTO_TEST_CASE(default_implementation) {
  ScriptInterfaceTest si_test;

  map<string, Variant> test_parameters;

  vector<int> vec_int{1, 4, 7, 10};
  vector<double> vec_double{0.1, 10.4, 11.9, 14.3};

  /* Some parameters to set, types should not
   * be relevant. */
  test_parameters["bool_opt"] = false;
  test_parameters["bool_req"] = true;
  test_parameters["integer"] = 5;
  test_parameters["vec_double"] = vec_double;
  test_parameters["vec_int"] = vec_int;

  /* Set values using the default implementation
   * set_parameters which in turn calls set_parameter
   * for every parameter. */
  si_test.set_parameters(test_parameters);

  /* Get values using the default implementation
   * set_parameter which in turn calls get_parameters
   * and extracts the required parameter. */
  BOOST_CHECK(boost::get<bool>(si_test.get_parameter("bool_opt")) == false);
  BOOST_CHECK(boost::get<bool>(si_test.get_parameter("bool_req")) == true);
  BOOST_CHECK(boost::get<int>(si_test.get_parameter("integer")) == 5);
  BOOST_CHECK(boost::get<vector<int>>(si_test.get_parameter("vec_int")) ==
              vec_int);
  BOOST_CHECK(boost::get<vector<double>>(si_test.get_parameter("vec_double")) ==
              vec_double);

  si_test.call_method("test_method", {});
}
