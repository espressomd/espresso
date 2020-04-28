/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef SCRIPT_INTERFACE_CYLINDER_WALL_HPP
#define SCRIPT_INTERFACE_CYLINDER_WALL_HPP

#include "Shape.hpp"
#include <shapes/Cylinder.hpp>

namespace ScriptInterface {
namespace Shapes {

class Cylinder : public Shape {
  using CoreShape = ::Shapes::Cylinder;
  std::shared_ptr<::Shapes::Cylinder> m_cylinder;

public:
  Cylinder() : m_cylinder(new ::Shapes::Cylinder()) {
    add_parameters(
        {{"radius", m_cylinder, &CoreShape::set_radius, &CoreShape::radius},
         {"length", m_cylinder, &CoreShape::set_length, &CoreShape::length},
         {"axis", m_cylinder, &CoreShape::set_axis, &CoreShape::axis},
         {"center", m_cylinder, &CoreShape::center},
         {"direction", m_cylinder, &CoreShape::direction},
         {"open", m_cylinder, &CoreShape::open}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_cylinder; }
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
