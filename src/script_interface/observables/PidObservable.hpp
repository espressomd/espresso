
/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "Observable.hpp"
#include "core/observables/PidObservable.hpp"

#include <memory>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace Observables {

/** Base class for script interfaces to particle-based observables
 *  @tparam CorePidObs Any core class derived from @ref
 *                     ::Observables::PidObservable "Observables::PidObservable"
 */
template <typename CorePidObs>
class PidObservable
    : public AutoParameters<PidObservable<CorePidObs>, Observable> {
public:
  static_assert(
      std::is_base_of<::Observables::PidObservable, CorePidObs>::value, "");

  PidObservable() {
    this->add_parameters({{"ids",
                           [this](Variant const &v) {
                             m_observable->ids() =
                                 get_value<std::vector<int>>(v);
                           },
                           [this]() { return m_observable->ids(); }}});
  }

  void construct(VariantMap const &params) override {
    m_observable =
        make_shared_from_args<CorePidObs, std::vector<int>>(params, "ids");
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
