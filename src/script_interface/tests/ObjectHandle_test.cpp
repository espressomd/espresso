/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE ObjectHandle test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/ObjectHandle.hpp"
#include "script_interface/ObjectState.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <utils/serialization/pack.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/variant.hpp>

#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

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
  static_assert(!std::is_copy_constructible_v<ObjectHandle>);
  static_assert(!std::is_copy_assignable_v<ObjectHandle>);
  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(default_constructible) {
  ObjectHandle handle;

  auto const param_name = std::string("unknown");
  handle.construct({{param_name, Variant{1}}});
  BOOST_CHECK(is_type<None>(handle.get_parameter(param_name)));
  BOOST_CHECK(is_type<None>(handle.call_method("foo", {})));

  // serialization should be empty
  auto const bytestring_obj = handle.serialize();
  auto const bytestring_ref = Utils::pack(ObjectState{});
  BOOST_CHECK_EQUAL(bytestring_obj, bytestring_ref);

  // serialization of an empty ObjectState should only contain the library
  // version and a string "serialization::archive", followed by a few integers
  auto const trim_null_terminator_right = [](std::string const &s) {
    return boost::trim_right_copy_if(s, [](char const c) { return c == '\0'; });
  };
  auto const bytestring_nul = Utils::pack(std::string{});
  auto const metadata_obj = trim_null_terminator_right(bytestring_obj);
  auto const metadata_ref = trim_null_terminator_right(bytestring_nul);
  BOOST_CHECK_EQUAL(metadata_obj, metadata_ref);
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

namespace boost {
namespace mpi {
class communicator {};
} // namespace mpi
} // namespace boost

namespace Testing {
/**
 * Logging mock for Context.
 */
struct LogContext : public Context {
  std::vector<std::pair<const ObjectHandle *, MockCall::Info>> call_log;

  void notify_call_method(const ObjectHandle *o, std::string const &n,
                          VariantMap const &p) override {
    call_log.emplace_back(o, MockCall::CallMethod{&n, &p});
  }
  void notify_set_parameter(const ObjectHandle *o, std::string const &n,
                            Variant const &v) override {
    call_log.emplace_back(o, MockCall::SetParameter{&n, &v});
  }

  std::shared_ptr<ObjectHandle> make_shared(std::string const &,
                                            const VariantMap &) override {
    auto it = std::make_shared<Testing::LogHandle>();
    set_context(it.get());

    return it;
  }
  std::shared_ptr<ObjectHandle>
  make_shared_local(std::string const &s, VariantMap const &v) override {
    return make_shared(s, v);
  }

  boost::string_ref name(const ObjectHandle *o) const override {
    return "Dummy";
  }

  bool is_head_node() const override { return true; }
  void parallel_try_catch(std::function<void()> const &) const override {}
  boost::mpi::communicator const &get_comm() const override { return m_comm; }

private:
  boost::mpi::communicator m_comm;
};
} // namespace Testing

/*
 * Check that Objecthandle::set_parameter does
 * notify the context.
 */
BOOST_AUTO_TEST_CASE(notify_set_parameter_) {
  using namespace Testing;
  auto log_ctx = std::make_shared<Testing::LogContext>();

  auto o = log_ctx->make_shared({}, {});

  std::string name;
  Variant value;

  o->set_parameter(name, value);

  auto const log_entry = log_ctx->call_log.at(0);
  BOOST_CHECK_EQUAL(log_entry.first, o.get());

  BOOST_CHECK((boost::get<MockCall::SetParameter>(log_entry.second) ==
               MockCall::SetParameter{&name, &value}));
}

/*
 * Check that Objecthandle::call_method does
 * notify the context.
 */
BOOST_AUTO_TEST_CASE(notify_call_method_) {
  using namespace Testing;
  auto log_ctx = std::make_shared<Testing::LogContext>();

  auto o = log_ctx->make_shared({}, {});

  std::string name;
  VariantMap params;
  o->call_method(name, params);

  auto const log_entry = log_ctx->call_log.at(0);
  BOOST_CHECK_EQUAL(log_entry.first, o.get());
  BOOST_CHECK((boost::get<MockCall::CallMethod>(log_entry.second) ==
               MockCall::CallMethod{&name, &params}));
}

/*
 * Check basic interface.
 */
BOOST_AUTO_TEST_CASE(interface_) {
  using namespace Testing;
  auto log_ctx = std::make_shared<Testing::LogContext>();
  auto o = log_ctx->make_shared({}, {});
  auto l = log_ctx->make_shared_local({}, {});
  BOOST_CHECK(log_ctx->is_head_node());
  BOOST_CHECK_EQUAL(log_ctx->name(o.get()), "Dummy");
  BOOST_CHECK_EQUAL(log_ctx->name(l.get()), "Dummy");
  log_ctx->parallel_try_catch([]() {});
  std::ignore = log_ctx->get_comm();
}
