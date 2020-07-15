/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#define BOOST_TEST_MODULE ObjectHandle test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/variant.hpp>

#include "script_interface/ScriptInterface.hpp"

using namespace ScriptInterface;

namespace Testing {
/* Types to keep track on calls on our mock handle.*/
namespace MockCall {
struct Construct {
  const VariantMap *params;

  bool operator==(Construct const &rhs) const { return params == rhs.params; }
};
struct SetParameter {
  const std::string *name;
  const Variant *value;

  bool operator==(SetParameter const &rhs) const {
    return std::tie(name, value) == std::tie(rhs.name, rhs.value);
  }
};
struct CallMethod {
  const std::string *name;
  const VariantMap *parameters;
  bool operator==(CallMethod const &rhs) const {
    return std::tie(name, parameters) == std::tie(rhs.name, rhs.parameters);
  }
};

using Info = boost::variant<Construct, SetParameter, CallMethod>;
} // namespace MockCall

/**
 * @brief Mock Objecthandle that logs all method calls.
 */
struct LogHandle : public ObjectHandle {
  std::vector<MockCall::Info> call_log;

  void do_construct(VariantMap const &params) override {
    call_log.emplace_back(MockCall::Construct{&params});
  }

  void do_set_parameter(const std::string &name,
                        const Variant &value) override {
    call_log.emplace_back(MockCall::SetParameter{&name, &value});
  }

  Variant do_call_method(const std::string &name,
                         const VariantMap &params) override {
    call_log.emplace_back(MockCall::CallMethod{&name, &params});

    return none;
  }
};
} // namespace Testing

BOOST_AUTO_TEST_CASE(non_copyable) {
  static_assert(!std::is_copy_constructible<ObjectHandle>::value, "");
  static_assert(!std::is_copy_assignable<ObjectHandle>::value, "");
}

/**
 * We check the default implementations of set_parameters
 * and get_parameter of ScriptInterface (this is the only
 * logic in the class).
 */
BOOST_AUTO_TEST_CASE(default_implementation) {
  Testing::LogHandle si_test;

  si_test.do_call_method("test_method", {});
}

/*
 * Check that the call to ObjectHandle::construct is
 * forwarded correctly to the implementation.
 */
BOOST_AUTO_TEST_CASE(do_construct_) {
  using namespace Testing;
  LogHandle log_handle;
  VariantMap test_params;

  log_handle.construct(test_params);
  BOOST_CHECK(boost::get<MockCall::Construct>(log_handle.call_log[0]) ==
              MockCall::Construct{&test_params});
}

/*
 * Check that the call to ObjectHandle::set_parameter is
 * forwarded correctly to the implementation.
 */
BOOST_AUTO_TEST_CASE(do_set_parameter_) {
  using namespace Testing;
  LogHandle log_handle;
  std::string name;
  Variant value;

  log_handle.set_parameter(name, value);
  BOOST_CHECK((boost::get<MockCall::SetParameter>(log_handle.call_log[0]) ==
               MockCall::SetParameter{&name, &value}));
}

/*
 * Check that the call to ObjectHandle::call_method is
 * forwarded correctly to the implementation.
 */
BOOST_AUTO_TEST_CASE(do_call_method_) {
  using namespace Testing;
  LogHandle log_handle;
  std::string name;
  VariantMap params;

  log_handle.call_method(name, params);
  BOOST_CHECK((boost::get<MockCall::CallMethod>(log_handle.call_log[0]) ==
               MockCall::CallMethod{&name, &params}));
}
