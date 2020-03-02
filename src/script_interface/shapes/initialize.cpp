/*
 * Copyright (C) 2015-2019 The ESPResSo project
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

#include "initialize.hpp"
#include "Cylinder.hpp"
#include "Ellipsoid.hpp"
#include "HollowConicalFrustum.hpp"
#include "NoWhere.hpp"
#include "Rhomboid.hpp"
#include "SimplePore.hpp"
#include "Slitpore.hpp"
#include "Sphere.hpp"
#include "SpheroCylinder.hpp"
#include "Stomatocyte.hpp"
#include "Torus.hpp"
#include "Union.hpp"
#include "Wall.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface {
namespace Shapes {
void initialize() {
  ScriptInterface::register_new<ScriptInterface::Shapes::HollowConicalFrustum>(
      "Shapes::HollowConicalFrustum");
  ScriptInterface::register_new<ScriptInterface::Shapes::Union>(
      "Shapes::Union");
  ScriptInterface::register_new<ScriptInterface::Shapes::NoWhere>(
      "Shapes::NoWhere");
  ScriptInterface::register_new<ScriptInterface::Shapes::Wall>("Shapes::Wall");
  ScriptInterface::register_new<ScriptInterface::Shapes::Ellipsoid>(
      "Shapes::Ellipsoid");
  ScriptInterface::register_new<ScriptInterface::Shapes::Sphere>(
      "Shapes::Sphere");
  ScriptInterface::register_new<ScriptInterface::Shapes::Cylinder>(
      "Shapes::Cylinder");
  ScriptInterface::register_new<ScriptInterface::Shapes::SpheroCylinder>(
      "Shapes::SpheroCylinder");
  ScriptInterface::register_new<ScriptInterface::Shapes::Rhomboid>(
      "Shapes::Rhomboid");
  ScriptInterface::register_new<ScriptInterface::Shapes::Slitpore>(
      "Shapes::Slitpore");
  ScriptInterface::register_new<ScriptInterface::Shapes::Stomatocyte>(
      "Shapes::Stomatocyte");
  ScriptInterface::register_new<ScriptInterface::Shapes::SimplePore>(
      "Shapes::SimplePore");
  ScriptInterface::register_new<ScriptInterface::Shapes::Torus>(
      "Shapes::Torus");
}
} /* namespace Shapes */
} /* namespace ScriptInterface */
