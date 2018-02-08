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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_AUTOUPDATEOBSERVABLES_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_AUTOUPDATEOBSERVABLES_HPP

#include "Observable.hpp"
#include "ScriptInterface.hpp"
#include "ScriptObjectRegistry.hpp"
#include "core/observables.hpp"

namespace ScriptInterface {
namespace Observables {

class AutoUpdateObservables : public ScriptObjectRegistry<Observable> {
  virtual void add_in_core(std::shared_ptr<Observable> obj_ptr) override {
    ::Observables::auto_update_observables.push_back(obj_ptr->observable());
  }
  virtual void remove_in_core(std::shared_ptr<Observable> obj_ptr) override {
    auto it = std::find(::Observables::auto_update_observables.begin(),
                        ::Observables::auto_update_observables.end(),
                        obj_ptr->observable());
    if (it != ::Observables::auto_update_observables.end()) {
      ::Observables::auto_update_observables.erase(it);
    }
  }
};
} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
