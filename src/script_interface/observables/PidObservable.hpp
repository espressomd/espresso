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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP

#include "ScriptInterface.hpp"

#include "Observable.hpp"
#include "core/observables/PidObservable.hpp"

#include <memory>
#include <type_traits>

namespace ScriptInterface {
namespace Observables {

template <typename CorePidObs> class PidObservable : public Observable {
public:
  static_assert(
      std::is_base_of<::Observables::PidObservable, CorePidObs>::value, "");
  
  PidObservable() : m_observable(std::make_shared<CorePidObs>()) {}

  VariantMap get_parameters() const override {
    return {{"ids", m_observable->ids()}};
  }

  ParameterMap valid_parameters() const override {
    return {{"ids", {ParameterType::INT_VECTOR, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if ("ids" == name) {
      m_observable->ids() = get_value<std::vector<int>>(value);
    }
  }

  virtual std::shared_ptr<::Observables::PidObservable> pid_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CorePidObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
