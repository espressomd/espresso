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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PARAMLESSOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PARAMLESSOBSERVABLE_HPP

#include "config.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "Observable.hpp"
#ifdef DPD
#include "core/observables/DPDStress.hpp"
#endif
#include "core/observables/LBFluidStress.hpp"
#include "core/observables/Observable.hpp"
#include "core/observables/StressTensor.hpp"

#include <memory>

namespace ScriptInterface {
namespace Observables {

/** Cython interface for parameter-free observables.
 *  Create new observables with @ref NEW_PARAMLESS_OBSERVABLE.
 *  @tparam CoreObs The core class exposed in Cython by this class.
 */
template <class CoreObs>
class ParamlessObservableInterface : public Observable {
public:
  ParamlessObservableInterface() : m_observable(new CoreObs()) {}

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

#define NEW_PARAMLESS_OBSERVABLE(name)                                         \
  using name = ParamlessObservableInterface<::Observables::name>;
NEW_PARAMLESS_OBSERVABLE(StressTensor)
NEW_PARAMLESS_OBSERVABLE(LBFluidStress)
#ifdef DPD
NEW_PARAMLESS_OBSERVABLE(DPDStress)
#endif

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
