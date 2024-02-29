/*
 * Copyright (C) 2022 The ESPResSo project
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
#define BOOST_TEST_MODULE Accumulators test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/variant.hpp>

#include "script_interface/GlobalContext.hpp"

#include "script_interface/accumulators/Correlator.hpp"
#include "script_interface/accumulators/MeanVarianceCalculator.hpp"
#include "script_interface/accumulators/TimeSeries.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/observables/ParamlessObservable.hpp"

#include "core/observables/Observable.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/errorhandling.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

static std::weak_ptr<boost::mpi::environment> mpi_env;

namespace Observables {
class MockObservable : public Observable {
public:
  std::vector<double>
  operator()(boost::mpi::communicator const &comm) const override {
    return {1., 2., 3., 4.};
  }
  std::vector<std::size_t> shape() const override { return {2u, 2u}; }
};
} // namespace Observables

namespace ScriptInterface {
namespace Observables {
NEW_PARAMLESS_OBSERVABLE(MockObservable)
} // namespace Observables
} // namespace ScriptInterface

using TestObs = ScriptInterface::Observables::MockObservable;
using TestObsPtr = std::shared_ptr<TestObs>;
using TimeSeries = ScriptInterface::Accumulators::TimeSeries;
using Correlator = ScriptInterface::Accumulators::Correlator;
using MeanVarianceCalculator =
    ScriptInterface::Accumulators::MeanVarianceCalculator;

using namespace ScriptInterface;

auto make_global_context(std::shared_ptr<Communication::MpiCallbacks> &cb) {
  Utils::Factory<ScriptInterface::ObjectHandle> factory;
  factory.register_new<TestObs>("TestObs");
  factory.register_new<TimeSeries>("TimeSeries");
  factory.register_new<Correlator>("Correlator");
  factory.register_new<MeanVarianceCalculator>("MeanVarianceCalculator");

  return std::make_shared<ScriptInterface::GlobalContext>(
      cb, std::make_shared<ScriptInterface::LocalContext>(factory, cb->comm()));
}

BOOST_AUTO_TEST_CASE(time_series) {
  boost::mpi::communicator world;
  auto cb =
      std::make_shared<Communication::MpiCallbacks>(world, ::mpi_env.lock());
  auto ctx = make_global_context(cb);

  if (world.rank() != 0) {
    cb->loop();
  } else {
    auto const obs = ctx->make_shared("TestObs", {});
    BOOST_REQUIRE(obs != nullptr);
    BOOST_CHECK_EQUAL(obs->context(), ctx.get());

    auto const acc_si =
        ctx->make_shared("TimeSeries", {{"obs", obs}, {"delta_N", 2}});
    BOOST_REQUIRE(acc_si != nullptr);
    auto const acc = std::dynamic_pointer_cast<TimeSeries>(acc_si);
    BOOST_REQUIRE(acc != nullptr);

    acc_si->call_method("update", VariantMap{});
    {
      BOOST_CHECK_EQUAL(get_value<int>(acc_si->get_parameter("delta_N")), 2);
      BOOST_CHECK_EQUAL(get_value<TestObsPtr>(acc_si->get_parameter("obs")),
                        obs);
    }
    {
      auto const shape_ref = std::vector<int>{{1, 2, 2}};
      // check non-const access
      auto const variant = acc_si->call_method("shape", VariantMap{});
      auto const shape = get_value<std::vector<int>>(variant);
      BOOST_TEST(shape == shape_ref, boost::test_tools::per_element());
      // check const access
      auto const shape_const = std::as_const(*acc).accumulator()->shape();
      BOOST_TEST(shape_const == shape_ref, boost::test_tools::per_element());
    }
    {
      auto const variant = acc_si->call_method("time_series", VariantMap{});
      auto const time_series = get_value<std::vector<Variant>>(variant);
      BOOST_REQUIRE_EQUAL(time_series.size(), 1u);
      auto const series = get_value<std::vector<double>>(time_series[0]);
      auto const series_ref = std::vector<double>{1., 2., 3., 4.};
      BOOST_TEST(series == series_ref, boost::test_tools::per_element());
    }
  }
}

BOOST_AUTO_TEST_CASE(correlator) {
  boost::mpi::communicator world;
  auto cb =
      std::make_shared<Communication::MpiCallbacks>(world, ::mpi_env.lock());
  auto ctx = make_global_context(cb);

  if (world.rank() != 0) {
    cb->loop();
  } else {
    auto const obs = ctx->make_shared("TestObs", {});
    BOOST_REQUIRE(obs != nullptr);

    auto const acc_si = ctx->make_shared(
        "Correlator",
        {{"obs1", obs},
         {"delta_N", 2},
         {"tau_lin", 4},
         {"tau_max", 2.},
         {"corr_operation", std::string("componentwise_product")}});
    BOOST_REQUIRE(acc_si != nullptr);
    auto acc = std::dynamic_pointer_cast<Correlator>(acc_si);
    BOOST_REQUIRE(acc != nullptr);

    acc_si->call_method("update", VariantMap{});
    acc_si->call_method("finalize", VariantMap{});
    {
      BOOST_CHECK_EQUAL(get_value<int>(acc_si->get_parameter("delta_N")), 2);
      BOOST_CHECK_EQUAL(get_value<int>(acc_si->get_parameter("tau_lin")), 4);
      BOOST_CHECK_EQUAL(get_value<double>(acc_si->get_parameter("tau_max")),
                        2.);
      BOOST_CHECK_EQUAL(
          get_value<std::string>(acc_si->get_parameter("corr_operation")),
          std::string("componentwise_product"));
      BOOST_CHECK_EQUAL(get_value<TestObsPtr>(acc_si->get_parameter("obs1")),
                        obs);
      BOOST_CHECK_EQUAL(get_value<TestObsPtr>(acc_si->get_parameter("obs2")),
                        obs);
    }
    {
      auto const shape_ref = std::vector<int>{{5, 2, 2}};
      // check non-const access
      auto const variant = acc_si->call_method("shape", VariantMap{});
      auto const shape = get_value<std::vector<int>>(variant);
      BOOST_TEST(shape == shape_ref, boost::test_tools::per_element());
      // check const access
      auto const shape_const = std::as_const(*acc).accumulator()->shape();
      BOOST_TEST(shape_const == shape_ref, boost::test_tools::per_element());
    }
    {
      auto const variant = acc_si->call_method("get_correlation", VariantMap{});
      auto const correlation = get_value<std::vector<double>>(variant);
      BOOST_REQUIRE_EQUAL(correlation.size(), 20u);
      auto corr_ref = std::vector<double>{1., 4., 9., 16.};
      corr_ref.resize(correlation.size());
      BOOST_TEST(correlation == corr_ref, boost::test_tools::per_element());
    }
    {
      auto const variant = acc_si->call_method("get_lag_times", VariantMap{});
      auto const lag_times = get_value<std::vector<double>>(variant);
      BOOST_REQUIRE_EQUAL(lag_times.size(), 5u);
      auto const lag_times_ref = std::vector<double>{0., -2., -4., -6., -8.};
      BOOST_TEST(lag_times == lag_times_ref, boost::test_tools::per_element());
    }
    {
      auto const variant =
          acc_si->call_method("get_samples_sizes", VariantMap{});
      auto const samples_n = get_value<std::vector<int>>(variant);
      BOOST_REQUIRE_EQUAL(samples_n.size(), 5u);
      auto const samples_n_ref = std::vector<int>{1, 0, 0, 0, 0};
      BOOST_TEST(samples_n == samples_n_ref, boost::test_tools::per_element());
    }
  }
}

BOOST_AUTO_TEST_CASE(mean_variance) {
  boost::mpi::communicator world;
  auto cb =
      std::make_shared<Communication::MpiCallbacks>(world, ::mpi_env.lock());
  auto ctx = make_global_context(cb);

  if (world.rank() != 0) {
    cb->loop();
  } else {
    auto const obs = ctx->make_shared("TestObs", {});

    auto const acc_si = ctx->make_shared("MeanVarianceCalculator",
                                         {{"obs", obs}, {"delta_N", 2}});
    BOOST_REQUIRE(acc_si != nullptr);
    auto acc = std::dynamic_pointer_cast<MeanVarianceCalculator>(acc_si);
    BOOST_REQUIRE(acc != nullptr);

    acc_si->call_method("update", VariantMap{});
    acc_si->call_method("update", VariantMap{});
    {
      BOOST_CHECK_EQUAL(get_value<int>(acc_si->get_parameter("delta_N")), 2);
      BOOST_CHECK_EQUAL(get_value<TestObsPtr>(acc_si->get_parameter("obs")),
                        obs);
    }
    {
      auto const shape_ref = std::vector<int>{{2, 2}};
      // check non-const access
      auto const variant = acc->call_method("shape", VariantMap{});
      auto const shape = get_value<std::vector<int>>(variant);
      BOOST_TEST(shape == shape_ref, boost::test_tools::per_element());
      // check const access
      auto const shape_const = std::as_const(*acc).accumulator()->shape();
      BOOST_TEST(shape_const == shape_ref, boost::test_tools::per_element());
    }
    {
      auto const variant = acc_si->call_method("mean", VariantMap{});
      auto const mean = get_value<std::vector<double>>(variant);
      BOOST_REQUIRE_EQUAL(mean.size(), 4u);
      auto const mean_ref = std::vector<double>{1., 2., 3., 4.};
      BOOST_TEST(mean == mean_ref, boost::test_tools::per_element());
    }
    {
      auto const variant = acc_si->call_method("variance", VariantMap{});
      auto const variance = get_value<std::vector<double>>(variant);
      BOOST_REQUIRE_EQUAL(variance.size(), 4u);
      auto const variance_ref = std::vector<double>{0., 0., 0., 0.};
      BOOST_TEST(variance == variance_ref, boost::test_tools::per_element());
    }
    {
      auto const variant = acc_si->call_method("std_error", VariantMap{});
      auto const stderror = get_value<std::vector<double>>(variant);
      BOOST_REQUIRE_EQUAL(stderror.size(), 4u);
      auto const stderror_ref = std::vector<double>{0., 0., 0., 0.};
      BOOST_TEST(stderror == stderror_ref, boost::test_tools::per_element());
    }
  }
}

int main(int argc, char **argv) {
  auto const mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  ::mpi_env = mpi_env;

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
