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

#ifndef SCRIPT_INTERFACE_CORRELATORS_CORRELATOR_HPP
#define SCRIPT_INTERFACE_CORRELATORS_CORRELATOR_HPP

#include "AccumulatorBase.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/accumulators/Correlator.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <string>
#include <utility>

namespace ScriptInterface {
namespace Accumulators {

class Correlator : public AccumulatorBase {
  using CoreCorr = ::Accumulators::Correlator;

public:
  Correlator() {
    /* Only args can be changed after construction. */
    add_parameters(
        {{"tau_lin", m_correlator, &CoreCorr::tau_lin},
         {"tau_max", m_correlator, &CoreCorr::tau_max},
         {"compress1", m_correlator, &CoreCorr::compress1},
         {"compress2", m_correlator, &CoreCorr::compress2},
         {"corr_operation", m_correlator, &CoreCorr::correlation_operation},
         {"args", m_correlator, &CoreCorr::set_correlation_args,
          &CoreCorr::correlation_args},
         {"obs1", std::as_const(m_obs1)},
         {"obs2", std::as_const(m_obs2)}});
  }

  void do_construct(VariantMap const &args) override {
    set_from_args(m_obs1, args, "obs1");
    if (args.count("obs2"))
      set_from_args(m_obs2, args, "obs2");
    else
      m_obs2 = m_obs1;

    auto const comp1 = get_value_or<std::string>(args, "compress1", "discard2");
    auto const comp2 = get_value_or<std::string>(args, "compress2", comp1);

    m_correlator = std::make_shared<CoreCorr>(
        get_value<int>(args, "tau_lin"), get_value<double>(args, "tau_max"),
        get_value<int>(args, "delta_N"), comp1, comp2,
        get_value<std::string>(args, "corr_operation"), m_obs1->observable(),
        m_obs2->observable(), get_value_or<Utils::Vector3d>(args, "args", {}));
  }

  std::shared_ptr<::Accumulators::Correlator> correlator() {
    return m_correlator;
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "update")
      correlator()->update();
    if (method == "finalize")
      correlator()->finalize();
    if (method == "get_correlation")
      return correlator()->get_correlation();
    if (method == "get_lag_times")
      return correlator()->get_lag_times();
    if (method == "get_samples_sizes")
      return correlator()->get_samples_sizes();

    return AccumulatorBase::call_method(method, parameters);
  }

  std::shared_ptr<::Accumulators::AccumulatorBase> accumulator() override {
    return m_correlator;
  }

  std::shared_ptr<const ::Accumulators::AccumulatorBase>
  accumulator() const override {
    return std::static_pointer_cast<::Accumulators::AccumulatorBase>(
        m_correlator);
  }

private:
  /* The actual correlator */
  std::shared_ptr<CoreCorr> m_correlator;

  std::shared_ptr<Observables::Observable> m_obs1;
  std::shared_ptr<Observables::Observable> m_obs2;

  std::string get_internal_state() const override {
    return m_correlator->get_internal_state();
  }

  void set_internal_state(std::string const &state) override {
    m_correlator->set_internal_state(state);
  }
};

} // namespace Accumulators
} /* namespace ScriptInterface */

#endif
