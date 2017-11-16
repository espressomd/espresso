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

#ifndef SCRIPT_INTERFACE_CORRELATORS_AUTOUPDATEACCUMULATORS_HPP
#define SCRIPT_INTERFACE_CORRELATORS_AUTOUPDATEACCUMULATORS_HPP

#include "Accumulator.hpp"
#include "ScriptInterface.hpp"
#include "ScriptObjectRegistry.hpp"
#include "core/accumulators.hpp"

namespace ScriptInterface {
namespace Accumulators {

class AutoUpdateAccumulators : public ScriptObjectRegistry<Correlator> {
  virtual void add_in_core(std::shared_ptr<Correlator> obj_ptr) {
    obj_ptr->accumulator()->start_auto_update();
    ::Accumulators::auto_update_accumulators.push_back(obj_ptr->accumulator());
  }
  virtual void remove_in_core(std::shared_ptr<Correlator> obj_ptr) {
    auto it = std::find(::Accumulators::auto_update_accumulators.begin(),
                        ::Accumulators::auto_update_accumulators.end(),
                        obj_ptr->accumulator());

    if (it != ::Accumulators::auto_update_accumulators.end()) {
      obj_ptr->accumulator()->stop_auto_update();
      ::Accumulators::auto_update_accumulators.erase(it);

    } else {
      throw "Could not find Accumulator to remove";
    };
  };

public:
  virtual const std::string name() const override {
    return "Accumulators::AutoUpdateAccumulators";
  };
};
} /* namespace Accumulators */
} /* namespace ScriptInterface */

#endif
