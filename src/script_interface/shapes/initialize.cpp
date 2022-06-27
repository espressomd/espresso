/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

#include "Cylinder.hpp"
#include "Ellipsoid.hpp"
#include "HollowConicalFrustum.hpp"
#include "NoWhere.hpp"
#include "Rhomboid.hpp"
#include "SimplePore.hpp"
#include "Slitpore.hpp"
#include "Sphere.hpp"
#include "SpheroCylinder.hpp"
#include "Torus.hpp"
#include "Union.hpp"
#include "Wall.hpp"

#include "initialize.hpp"

namespace ScriptInterface {
namespace Shapes {
void initialize(Utils::Factory<ObjectHandle> *f) {
  f->register_new<HollowConicalFrustum>("Shapes::HollowConicalFrustum");
  f->register_new<Union>("Shapes::Union");
  f->register_new<NoWhere>("Shapes::NoWhere");
  f->register_new<Wall>("Shapes::Wall");
  f->register_new<Ellipsoid>("Shapes::Ellipsoid");
  f->register_new<Sphere>("Shapes::Sphere");
  f->register_new<Cylinder>("Shapes::Cylinder");
  f->register_new<SpheroCylinder>("Shapes::SpheroCylinder");
  f->register_new<Rhomboid>("Shapes::Rhomboid");
  f->register_new<Slitpore>("Shapes::Slitpore");
  f->register_new<SimplePore>("Shapes::SimplePore");
  f->register_new<Torus>("Shapes::Torus");
}
} /* namespace Shapes */
} /* namespace ScriptInterface */
