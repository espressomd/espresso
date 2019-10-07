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
#include "Torus.hpp"
#include "Union.hpp"
#include "Wall.hpp"

#include "script_interface/ClassName.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <utils/tuple.hpp>

#include <tuple>

namespace ScriptInterface {
namespace Shapes {
constexpr auto class_names() {
  return std::make_tuple(
      ClassName<HollowConicalFrustum>{"Shapes::HollowConicalFrustum"},
      ClassName<Union>{"Shapes::Union"}, ClassName<NoWhere>{"Shapes::NoWhere"},
      ClassName<Wall>{"Shapes::Wall"},
      ClassName<Ellipsoid>{"Shapes::Ellipsoid"},
      ClassName<Sphere>{"Shapes::Sphere"},
      ClassName<Cylinder>{"Shapes::Cylinder"},
      ClassName<SpheroCylinder>{"Shapes::SpheroCylinder"},
      ClassName<Rhomboid>{"Shapes::Rhomboid"},
      ClassName<Slitpore>{"Shapes::Slitpore"},
      ClassName<SimplePore>{"Shapes::SimplePore"},
      ClassName<Torus>{"Shapes::Torus"});
}

void initialize(ObjectManager *f) {
  Utils::for_each(
      [f](auto name) {
        f->register_new<typename decltype(name)::class_type>(name.name);
      },
      class_names());
}
} /* namespace Shapes */
} /* namespace ScriptInterface */
