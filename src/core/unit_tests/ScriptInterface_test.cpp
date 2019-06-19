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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ScriptInterface test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpi.hpp>

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/ObjectManager.hpp"

using std::map;
using std::string;
using std::vector;

using namespace ScriptInterface;

/**
 * @brief Mock to test ScriptInterface.
 */
struct ScriptInterfaceTest : public ScriptInterface::ObjectHandle {
  Variant get_parameter(const std::string &name) const override {
    if(name == "test_parameter") {
      return test_parameter;
    } else {
      return none;
    }
  }

  /* Not needed for testing */
  Utils::Span<const boost::string_ref> valid_parameters() const override {
    static std::array<const boost::string_ref, 1> params = {"test_parameter"};
    return params;
  }

  void do_set_parameter(const string &name, const Variant &value) override {
    if(name == "test_parameter") {
      test_parameter = boost::get<int>(value);
    }
  }

  Variant do_call_method(const std::string &name,
                         const VariantMap &params) override {
    if (name == "test_method") {
      test_method_called = true;
    }

    return none;
  }

  int test_parameter = -1;
  bool test_method_called = false;
};

BOOST_AUTO_TEST_CASE(non_copyable) {
  static_assert(!std::is_copy_constructible<ObjectHandle>::value, "");
  static_assert(!std::is_copy_assignable<ObjectHandle>::value, "");
}

/*
 * We check the default implementations of set_parameters
 * and get_parameter of ScriptInterface (this is the only
 * logic in the class).
 */
BOOST_AUTO_TEST_CASE(default_implementation) {
  ScriptInterfaceTest si_test;

  si_test.do_call_method("test_method", {});
}

BOOST_AUTO_TEST_CASE(set_parameter_test) {
  boost::mpi::communicator world;
  Communication::MpiCallbacks cb{world};
  ObjectManager om(&cb);

  if(world.rank() == 0) {
  } else {
    cb.loop();
  }

}

int main(int argc, char **argv) {
boost::mpi::environment mpi_env(argc, argv);

register_new<ScriptInterfaceTest>("TestClass");

return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
