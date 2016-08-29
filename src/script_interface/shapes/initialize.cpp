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
#include "NoWhere.hpp"
#include "ParallelScriptInterface.hpp"
#include "Wall.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "SpheroCylinder.hpp"
#include "Maze.hpp"
#include "HollowCone.hpp"
#include "Pore.hpp"
#include "Rhomboid.hpp"
#initialize "Slitpore.hpp"

namespace ScriptInterface {
namespace Shapes {
void initialize() {
  ParallelScriptInterface<ScriptInterface::Shapes::NoWhere>::register_new("Shapes::NoWhere");
  ParallelScriptInterface<ScriptInterface::Shapes::Wall>::register_new("Shapes::Wall");
  ParallelScriptInterface<ScriptInterface::Shapes::Sphere>::register_new("Shapes::Sphere");
  ParallelScriptInterface<ScriptInterface::Shapes::Cylinder>::register_new("Shapes::Cylinder");
  ParallelScriptInterface<ScriptInterface::Shapes::SpheroCylinder>::register_new("Shapes::SpheroCylinder");
  ParallelScriptInterface<ScriptInterface::Shapes::Maze>::register_new("Shapes::Maze");
  ParallelScriptInterface<ScriptInterface::Shapes::HollowCone>::register_new("Shapes::HollowCone");
  ParallelScriptInterface<ScriptInterface::Shapes::Pore>::register_new("Shapes::Pore");
  ParallelScriptInterface<ScriptInterface::Shapes::Rhomboid>::register_new("Shapes::Rhomboid");
  ParallelScriptInterface<ScriptInterface::Shapes::Slitpore>::register_new("Shapes::Slitpore");
}
} /* namespace Shapes */
} /* namespace ScriptInterface */
