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

#include <memory>

namespace ScriptInterface {
namespace Correlators {

class Correlator : public AutoParameters {
public:
  Correlator() : m_correlator(std::make_shared<::Correlators::Correlator>()) {
    add_parameters({{"tau_lin", m_correlator->tau_lin},
                    {"tau_max", m_correlator->tau_max},
                    {"dt", m_correlator->dt},
                    {"compress1", m_correlator->compressA_name},
                    {"compress2", m_correlator->compressB_name},
                    {"corr_operation", m_correlator->corr_operation_name},
                    {"args", m_correlator->correlation_args},
                    {"obs1", m_obs1},
                    {"obs2", m_obs2}});
  }

  std::shared_ptr<::Correlators::Correlator> correlator() {
    return m_correlator;
  }

  void check_if_initialized() {
    if (!m_correlator->initialized)
      throw std::runtime_error("The correlator has not yet been initialied.");
  }

  virtual Variant call_method(std::string const &method,
                              VariantMap const &parameters) override {
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
