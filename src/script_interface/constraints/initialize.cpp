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
#include "initialize.hpp"
#include "ScriptInterface.hpp"

#include "Constraints.hpp"

#include "ShapeBasedConstraint.hpp"
#include "HomogeneousMagneticField.hpp"

namespace ScriptInterface {
namespace Constraints {

void initialize() {
  ScriptInterface::register_new<ScriptInterface::Constraints::Constraints>(
      "Constraints::Constraints");

  ScriptInterface::register_new<ScriptInterface::Constraints::ShapeBasedConstraint>(
      "Constraints::ShapeBasedConstraint");

  ScriptInterface::register_new<ScriptInterface::Constraints::HomogeneousMagneticField>(
      "Constraints::HomogeneousMagneticField");
}
} /* namespace Constraints */
} /* namespace ScriptInterface */
