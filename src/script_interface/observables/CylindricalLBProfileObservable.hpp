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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP

#include "auto_parameters/AutoParameters.hpp"

#include <memory>

#include "CylindricalProfileObservable.hpp"
#include "LBObservable.hpp"
#include "Observable.hpp"
#include "core/observables/CylindricalLBProfileObservable.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename CoreCylLBObs>
class CylindricalLBProfileObservable
    : public CylindricalProfileObservable<CoreCylLBObs>,
      public LBObservable<CoreCylLBObs> {
public:
  static_assert(std::is_base_of<::Observables::CylindricalLBProfileObservable,
                                CoreCylLBObs>::value,
                "");
  CylindricalLBProfileObservable()
      : m_observable(std::make_shared<CoreCylLBObs>()) {}

  virtual Variant call_method(std::string const &method,
                              VariantMap const &parameters) override {
    if (method == "calculate") {
      return cylindrical_profile_observable()->operator()(partCfg());
    }
    return {};
  }

  virtual std::shared_ptr<::Observables::Observable>
  observable() const override {
    return m_observable;
  }

  virtual std::shared_ptr<::Observables::CylindricalLBProfileObservable>
  cylindrical_profile_observable() const {
    return m_observable;
  }

private:
  std::shared_ptr<CoreCylLBObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
