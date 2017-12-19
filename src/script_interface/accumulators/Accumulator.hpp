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

#ifndef SCRIPT_INTERFACE_ACCUMULATORS_ACCUMULATOR_HPP
#define SCRIPT_INTERFACE_ACCUMULATORS_ACCUMULATOR_HPP

#include "ScriptInterface.hpp"
#include "core/accumulators/Accumulator.hpp"
#include "core/utils/Factory.hpp"
#include "observables/Observable.hpp"

#include <memory>

namespace ScriptInterface {
namespace Accumulators {

class Accumulator : public AutoParameters {
public:
  Accumulator()
      : m_accumulator(std::make_shared<::Accumulators::Accumulator>())
  {
    add_parameters({{"obs",
                     [this](Variant const &value) {
                       auto obs_ptr = get_value<std::shared_ptr<Observables::Observable>>(value);
                       // We are expecting a ScriptInterface::Observables::Observable here,
                       // throw if not. That means the assigned object had the wrong type.
                       if (obs_ptr) {
                         m_obs = obs_ptr;
                         m_accumulator->m_obs = obs_ptr->observable();
                         m_accumulator->initialize();
                       }
                     },
                     [this]() {
                       return m_obs ? m_obs->id() : ObjectId();
                     }}});
  }

  const std::string name() const override { return "Accumulators::Accumulator"; }

  std::shared_ptr<::Accumulators::Accumulator> accumulator() {
    return m_accumulator;
  }

  void check_if_initialized() {
    if (!m_accumulator->m_initialized)
      throw std::runtime_error("The accumulator has not yet been initialied.");
  }
  virtual Variant call_method(std::string const &method,
                              VariantMap const &parameters) override {
    check_if_initialized();
    if (method == "update") {
      if (m_accumulator->m_autoupdate) {
        throw std::runtime_error(
            "auto_update is enabled for the accumulator. Cannot update manually");
      }
      return m_accumulator->update();
    }
    if (method == "auto_update") 
      return m_accumulator->m_autoupdate;
    if (method == "get_mean")
      return m_accumulator->get_mean();
    if (method == "get_variance")
      return m_accumulator->get_variance();
    return {};
  }

private:
  /* The actual accumulator */
  std::shared_ptr<::Accumulators::Accumulator> m_accumulator;
  std::shared_ptr<Observables::Observable> m_obs;
};

} /* namespace Accumulator */
} /* namespace ScriptInterface */

#endif
