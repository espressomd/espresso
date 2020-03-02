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

#ifndef SCRIPT_INTERFACE_SHAPES_SPHEROCYLINDER_HPP
#define SCRIPT_INTERFACE_SHAPES_SPHEROCYLINDER_HPP

#include "Shape.hpp"
#include <shapes/SpheroCylinder.hpp>

namespace ScriptInterface {
namespace Shapes {

class SpheroCylinder : public Shape {
  using CoreShape = ::Shapes::SpheroCylinder;
  std::shared_ptr<::Shapes::SpheroCylinder> m_spherocylinder;

public:
  SpheroCylinder() : m_spherocylinder(new ::Shapes::SpheroCylinder()) {
    add_parameters(
        {{"radius", m_spherocylinder, &CoreShape::set_radius,
          &CoreShape::radius},
         {"length", m_spherocylinder, &CoreShape::set_length,
          &CoreShape::length},
         {"axis", m_spherocylinder, &CoreShape::set_axis, &CoreShape::axis},
         {"direction", m_spherocylinder, &CoreShape::direction},
         {"center", m_spherocylinder, &CoreShape::center}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_spherocylinder;
  }
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
