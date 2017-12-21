/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_CORRELATORS_CORRELATOR_HPP
#define SCRIPT_INTERFACE_CORRELATORS_CORRELATOR_HPP

#include "ScriptInterface.hpp"
#include "auto_parameters/AutoParameters.hpp"

#include "core/correlators/Correlator.hpp"
#include "observables/Observable.hpp"

#include "core/utils/as_const.hpp"

#include <memory>

namespace ScriptInterface {
namespace Correlators {

class Correlator : public AutoParameters {
  using CoreCorr = ::Correlators::Correlator;

public:
  Correlator() {
    using Utils::as_const;
    /* Only args can be changed after construction. */
    add_parameters(
        {{"tau_lin", m_correlator, &CoreCorr::tau_lin},
         {"tau_max", m_correlator, &CoreCorr::tau_max},
         {"dt", m_correlator, &CoreCorr::dt},
         {"compress1", m_correlator, &CoreCorr::compress1},
         {"compress2", m_correlator, &CoreCorr::compress2},
         {"corr_operation", m_correlator, &CoreCorr::correlation_operation},
         {"args", m_correlator, &CoreCorr::set_correlation_args,
          &CoreCorr::correlation_args},
         {"dim_corr", m_correlator, &CoreCorr::dim_corr},
         {"obs1", as_const(m_obs1)},
         {"obs2", as_const(m_obs2)},
         {"n_result", m_correlator, &CoreCorr::n_result}});
  }

  void construct(VariantMap const &args) override {
    set_from_args(m_obs1, args, "obs1");
    if (args.count("obs2"))
      set_from_args(m_obs2, args, "obs2");
    else
      m_obs2 = m_obs1;

    m_correlator = std::make_shared<CoreCorr>(
        get_value<int>(args, "tau_lin"), get_value<double>(args, "tau_max"),
        get_value<double>(args, "dt"),
        /* These two are optional */
        get_value_or<std::string>(args, "compess1", ""),
        get_value_or<std::string>(args, "compess2", ""),
        get_value<std::string>(args, "corr_operation"), m_obs1->observable(),
        m_obs2->observable());
  }

  std::shared_ptr<::Correlators::Correlator> correlator() {
    return m_correlator;
  }

  virtual Variant call_method(std::string const &method,
                              VariantMap const &parameters) override {
    if (method == "update") {
      if (m_correlator->autoupdate) {
        throw std::runtime_error(
            "auto_update is enable for the correlator. Cannot update manually");
      }
      if (m_correlator->get_data()) {
        throw std::runtime_error("Correlator update failed");
      }
    }
    if (method == "auto_update") {
      return m_correlator->autoupdate;
    }
    if (method == "finalize")
      m_correlator->finalize();
    if (method == "get_correlation") {
      return m_correlator->get_correlation();
    }

    return {};
  }

  Variant get_state() const override {
    std::vector<Variant> state(2);
    state[0] = ScriptInterfaceBase::get_state();
    state[1] = m_correlator->get_internal_state();

    return state;
  }

private:
  void set_state(Variant const &state) override {
    auto const &state_vec = boost::get<std::vector<Variant>>(state);

    ScriptInterfaceBase::set_state(state_vec.at(0));
    m_correlator->set_internal_state(boost::get<std::string>(state_vec.at(1)));
  }

  /* The actual correlator */
  std::shared_ptr<CoreCorr> m_correlator;

  std::shared_ptr<Observables::Observable> m_obs1;
  std::shared_ptr<Observables::Observable> m_obs2;
};

} /* namespace Correlators */
} /* namespace ScriptInterface */

#endif
