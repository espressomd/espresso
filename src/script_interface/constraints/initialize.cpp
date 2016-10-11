/*
  Copyright (C) 2015,2016 The ESPResSo project

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
#include "config.hpp" 

#include "initialize.hpp"
#include "ParallelScriptInterface.hpp"
#include "utils/Factory.hpp"

#include "Constraint.hpp"
#include "Constraints.hpp"

namespace ScriptInterface {
namespace Constraints {

void initialize() {
#ifdef CONSTRAINTS
  ParallelScriptInterface<ScriptInterface::Constraints::Constraints>::
    register_new("Constraints::Constraints");

  ParallelScriptInterface<ScriptInterface::Constraints::Constraint>::
    register_new("Constraints::Constraint");
#endif
}
} /* namespace Constraints */
} /* namespace ScriptInterface */

