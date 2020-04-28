/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <type_traits>

#define BOOST_TEST_MODULE ScriptInterface test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/ScriptInterface.hpp"

using std::map;
using std::string;
using std::vector;

using namespace ScriptInterface;

namespace Testing {

/* Mock to test ScriptInterface. */
struct ScriptInterfaceTest : public ScriptInterface::ScriptInterfaceBase {
  VariantMap get_parameters() const override {
    VariantMap ret;

    ret["bool_opt"] = bool_opt;
    ret["integer"] = integer;
    ret["bool_req"] = bool_req;
    ret["vec_double"] = vec_double;
    ret["vec_int"] = vec_int;

    return ret;
  }

  /* Not needed for testing */
  Utils::Span<const boost::string_ref> valid_parameters() const override {
    return {};
  }

  void set_parameter(const string &name, const Variant &value) override {
    if (name == "bool_opt") {
      bool_opt =
          get_value<std::remove_reference<decltype(bool_opt)>::type>(value);
    };
    if (name == "integer") {
      integer =
          get_value<std::remove_reference<decltype(integer)>::type>(value);
    };
    if (name == "bool_req") {
      bool_req =
          get_value<std::remove_reference<decltype(bool_req)>::type>(value);
    };
    if (name == "vec_double") {
      vec_double =
          get_value<std::remove_reference<decltype(vec_double)>::type>(value);
    };
    if (name == "vec_int") {
      vec_int =
          get_value<std::remove_reference<decltype(vec_int)>::type>(value);
    };
  }

  Variant call_method(const std::string &name,
                      const VariantMap &params) override {
    if (name == "test_method") {
    }

    return true;
  }

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

/*
 * We check the default implementations of set_parameters
 * and get_parameter of ScriptInterface (this is the only
 * logic in the class).
 */
BOOST_AUTO_TEST_CASE(default_implementation) {
  ScriptInterfaceTest si_test;

  si_test.call_method("test_method", {});
}
