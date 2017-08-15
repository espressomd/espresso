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
#include "core/correlators/Correlator.hpp"
#include "core/utils/Factory.hpp"
#include "observables/Observable.hpp"

#include <memory>

namespace ScriptInterface {
namespace Correlators {

class Correlator : public ScriptInterfaceBase {
public:
  Correlator()
      : m_correlator(new ::Correlators::Correlator()), m_obs1(nullptr),
        m_obs2(nullptr) {}

  const std::string name() const override { return "Correlators::Correlator"; }

  VariantMap get_parameters() const override {
    return {{"tau_lin", (int)m_correlator->tau_lin},
            {"tau_max", (int)m_correlator->tau_max},
            {"dt", m_correlator->dt},
            {"compress1", m_correlator->compressA_name},
            {"compress2", m_correlator->compressB_name},
            {"corr_operation", m_correlator->corr_operation_name},
            {"args", m_correlator->correlation_args},
            {"obs1", (m_obs1 != nullptr) ? m_obs1->id() : ObjectId()},
            {"obs2", (m_obs2 != nullptr) ? m_obs2->id() : ObjectId()}};
  }

  ParameterMap valid_parameters() const override {
    return {{"tau_lin", {ParameterType::INT, true}},
            {"tau_max", {ParameterType::DOUBLE, true}},
            {"dt", {ParameterType::DOUBLE, true}},
            {"obs1", {ParameterType::OBJECTID, true}},
            {"obs2", {ParameterType::OBJECTID, true}},
            {"compress1", {ParameterType::STRING, true}},
            {"compress2", {ParameterType::STRING, true}},
            {"args", {ParameterType::DOUBLE_VECTOR, true}},
            {"corr_operation", {ParameterType::STRING, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if (m_correlator->initialized) {
      throw std::runtime_error(
          "Correlator cannot be changed after initial setup");
    }
    if ((name == "obs1") || (name == "obs2")) {
      auto obs_ptr = get_value<std::shared_ptr<Observables::Observable>>(value);

      /* We are expecting a ScriptInterface::Observables::Observable here,
         throw if not. That means the assigned object had the wrong type. */
      if (obs_ptr != nullptr) {
        if (name == "obs1") {
          m_correlator->A_obs = obs_ptr->observable();
          m_obs1 = obs_ptr;
        }
        if (name == "obs2") {
          m_correlator->B_obs = obs_ptr->observable();
          m_obs2 = obs_ptr;
        }
      }
    }

    SET_PARAMETER_HELPER("tau_lin", m_correlator->tau_lin);
    SET_PARAMETER_HELPER("tau_max", m_correlator->tau_max);
    SET_PARAMETER_HELPER("dt", m_correlator->dt);
    SET_PARAMETER_HELPER("args", m_correlator->correlation_args);
    SET_PARAMETER_HELPER("corr_operation", m_correlator->corr_operation_name);
    SET_PARAMETER_HELPER("compress1", m_correlator->compressA_name);
    SET_PARAMETER_HELPER("compress2", m_correlator->compressB_name);
  }
  virtual void set_parameters(const VariantMap &p) override {
    ScriptInterfaceBase::set_parameters(p);
    m_correlator->initialize();
  };

  std::shared_ptr<::Correlators::Correlator> correlator() {
    return m_correlator;
  }
  void check_if_initialized() {
    if (!m_correlator->initialized)
      throw std::runtime_error("The correlator has not yet been initialied.");
  }
  virtual Variant call_method(std::string const &method,
                              VariantMap const &parameters) {
    check_if_initialized();
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
    if (method == "n_results")
      return m_correlator->n_result;
    if (method == "dim_corr")
      return m_correlator->dim_corr;
    if (method == "hierarchy_depth")
      return m_correlator->hierarchy_depth;

    return {};
  }

private:
  /* The actual correlator */
  std::shared_ptr<::Correlators::Correlator> m_correlator;

  std::shared_ptr<Observables::Observable> m_obs1;
  std::shared_ptr<Observables::Observable> m_obs2;
};

} /* namespace Correlators */
} /* namespace ScriptInterface */

#endif
