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

#ifndef SCRIPT_INTERFACE_CORRELATORS_AUTOUPDATECORRELATORS_HPP
#define SCRIPT_INTERFACE_CORRELATORS_AUTOUPDATECORRELATORS_HPP

#include "Correlator.hpp"
#include "ScriptInterface.hpp"
#include "ScriptObjectRegistry.hpp"
#include "core/correlators.hpp"

namespace ScriptInterface {
namespace Correlators {

class AutoUpdateCorrelators : public ScriptObjectRegistry<Correlator> {
  virtual void add_in_core(std::shared_ptr<Correlator> obj_ptr) {
    obj_ptr->correlator()->start_auto_update();
    ::Correlators::auto_update_correlators.push_back(obj_ptr->correlator());
  }
  virtual void remove_in_core(std::shared_ptr<Correlator> obj_ptr) {
    auto it = std::find(::Correlators::auto_update_correlators.begin(),
                        ::Correlators::auto_update_correlators.end(),
                        obj_ptr->correlator());

    if (it != ::Correlators::auto_update_correlators.end()) {
      obj_ptr->correlator()->stop_auto_update();
      ::Correlators::auto_update_correlators.erase(it);

    } else {
      throw "Could not find Correlator to remove";
    };
  };

public:
  virtual const std::string name() const override {
    return "Correlators::AutoUpdateCorrelators";
  };
};
} /* namespace Correlators */
} /* namespace ScriptInterface */

#endif
