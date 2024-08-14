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

#pragma once

#include "AccumulatorBase.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/accumulators/MeanVarianceCalculator.hpp"
#include "core/system/System.hpp"

#include <memory>
#include <string>
#include <utility>

namespace ScriptInterface {
namespace Accumulators {

class MeanVarianceCalculator : public AccumulatorBase {
public:
  MeanVarianceCalculator() { add_parameters({{"obs", std::as_const(m_obs)}}); }

  void do_construct(VariantMap const &params) override {
    set_from_args(m_obs, params, "obs");

    if (m_obs)
      m_accumulator = std::make_shared<::Accumulators::MeanVarianceCalculator>(
          get_core_system_pointer(params),
          get_value_or<int>(params, "delta_N", 1), m_obs->observable());
  }

  std::shared_ptr<::Accumulators::MeanVarianceCalculator>
  mean_variance_calculator() {
    return m_accumulator;
  }

  std::shared_ptr<const ::Accumulators::MeanVarianceCalculator>
  mean_variance_calculator() const {
    return m_accumulator;
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "update") {
      ObjectHandle::context()->parallel_try_catch(
          [&]() { mean_variance_calculator()->update(context()->get_comm()); });
      return {};
    }
    if (method == "mean")
      return mean_variance_calculator()->mean();
    if (method == "variance")
      return mean_variance_calculator()->variance();
    if (method == "std_error")
      return mean_variance_calculator()->std_error();
    return AccumulatorBase::do_call_method(method, parameters);
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
  std::shared_ptr<::Accumulators::MeanVarianceCalculator> m_accumulator;
  std::shared_ptr<Observables::Observable> m_obs;
};

} // namespace Accumulators
} // namespace ScriptInterface
