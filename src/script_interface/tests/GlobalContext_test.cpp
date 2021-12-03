/*
 * Copyright (C) 2017-2020 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ScriptInterface::GlobalContext test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/mpi.hpp>

#include "script_interface/GlobalContext.hpp"

#include <algorithm>
#include <memory>
#include <string>

namespace si = ScriptInterface;

struct Dummy : si::ObjectHandle {
  si::VariantMap params;

  si::Variant get_parameter(std::string const &name) const override {
    return params.at(name);
  }

  void do_set_parameter(std::string const &name,
                        si::Variant const &val) override {
    params[name] = val;
  }

  Utils::Span<const boost::string_ref> valid_parameters() const override {
    static const boost::string_ref parameter_names[] = {"id", "object_param"};

    return Utils::make_const_span(parameter_names,
                                  std::min(params.size(), 2lu));
  }
};

auto make_global_context(Communication::MpiCallbacks &cb) {
  Utils::Factory<si::ObjectHandle> factory;
  factory.register_new<Dummy>("Dummy");

  return std::make_shared<si::GlobalContext>(
      cb, std::make_shared<si::LocalContext>(factory, 0));
}

BOOST_AUTO_TEST_CASE(GlobalContext_make_shared) {
  boost::mpi::communicator world;
  Communication::MpiCallbacks cb{world};
  auto ctx = make_global_context(cb);

  if (world.rank() == 0) {
    auto res = ctx->make_shared("Dummy", {});
    BOOST_REQUIRE(res != nullptr);
    BOOST_CHECK_EQUAL(res->context(), ctx.get());
    BOOST_CHECK_EQUAL(ctx->name(res.get()), "Dummy");
  } else {
    cb.loop();
  }
}

BOOST_AUTO_TEST_CASE(GlobalContext_serialization) {
  boost::mpi::communicator world;
  Communication::MpiCallbacks cb{world};
  auto ctx = make_global_context(cb);

  if (world.rank() == 0) {
    auto const serialized = [&]() {
      auto d1 = ctx->make_shared("Dummy", {});
      auto d2 = ctx->make_shared("Dummy", {});
      auto d3 = ctx->make_shared("Dummy", {});

      d1->set_parameter("object_param", d2);
      d1->set_parameter("id", 1);
      d2->set_parameter("object_param", d3);
      d2->set_parameter("id", 2);
      d3->set_parameter("id", 3);

      return d1->serialize();
    }();

    auto d1 = si::ObjectHandle::deserialize(serialized, *ctx);
    BOOST_REQUIRE(d1);
    BOOST_CHECK_EQUAL(boost::get<int>(d1->get_parameter("id")), 1);
    auto d2 = boost::get<si::ObjectRef>(d1->get_parameter("object_param"));
    BOOST_REQUIRE(d2);
    BOOST_CHECK_EQUAL(boost::get<int>(d2->get_parameter("id")), 2);
    auto d3 = boost::get<si::ObjectRef>(d2->get_parameter("object_param"));
    BOOST_REQUIRE(d3);
    BOOST_CHECK_EQUAL(boost::get<int>(d3->get_parameter("id")), 3);
  } else {
    cb.loop();
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
