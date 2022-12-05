/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_ACCUMULATORS_TIME_SERIES_HPP
#define SCRIPT_INTERFACE_ACCUMULATORS_TIME_SERIES_HPP

#include "AccumulatorBase.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/accumulators/TimeSeries.hpp"

#include <boost/range/algorithm/transform.hpp>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Accumulators {

class TimeSeries : public AccumulatorBase {
public:
  TimeSeries() { add_parameters({{"obs", std::as_const(m_obs)}}); }

  void do_construct(VariantMap const &params) override {
    set_from_args(m_obs, params, "obs");

    if (m_obs)
      m_accumulator = std::make_shared<::Accumulators::TimeSeries>(
          m_obs->observable(), get_value_or<int>(params, "delta_N", 1));
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "update") {
      m_accumulator->update();
    }
    if (method == "time_series") {
      auto const &series = m_accumulator->time_series();
      std::vector<Variant> ret(series.size());

      boost::transform(
          series, ret.begin(),
          [](std::vector<double> const &sample) { return sample; });

      return ret;
    }
    if (method == "clear") {
      m_accumulator->clear();
    }

    return AccumulatorBase::call_method(method, parameters);
  }

  std::shared_ptr<::Accumulators::AccumulatorBase> accumulator() override {
    return m_accumulator;
  }

  std::shared_ptr<const ::Accumulators::AccumulatorBase>
  accumulator() const override {
    return std::static_pointer_cast<::Accumulators::AccumulatorBase>(
        m_accumulator);
  }

private:
  /* The actual accumulator */
  std::shared_ptr<::Accumulators::TimeSeries> m_accumulator;
  std::shared_ptr<Observables::Observable> m_obs;

  std::string get_internal_state() const override {
    return m_accumulator->get_internal_state();
  }

  void set_internal_state(std::string const &state) override {
    m_accumulator->set_internal_state(state);
  }
};

} // namespace Accumulators
} /* namespace ScriptInterface */

#endif // SCRIPT_INTERFACE_ACCUMULATORS_TIMESERIES_HPP
