/*
 * Copyright (C) 2015-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "BindAtPointOfCollision.hpp"
#include "BindCenters.hpp"
#include "CollisionDetection.hpp"
#include "GlueToSurface.hpp"
#include "Off.hpp"
#include "initialize.hpp"

namespace ScriptInterface::CollisionDetection {
void initialize(Utils::Factory<ObjectHandle> *om) {
#ifdef ESPRESSO_COLLISION_DETECTION
  om->register_new<CollisionDetection>(
      "CollisionDetection::CollisionDetection");
  om->register_new<Off>("CollisionDetection::Off");
  om->register_new<BindCenters>("CollisionDetection::BindCenters");
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  om->register_new<BindAtPointOfCollision>(
      "CollisionDetection::BindAtPointOfCollision");
  om->register_new<GlueToSurface>("CollisionDetection::GlueToSurface");
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#endif // ESPRESSO_COLLISION_DETECTION
}
} // namespace ScriptInterface::CollisionDetection
