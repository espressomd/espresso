/*
  Copyright (C) 2017 The ESPResSo project

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

#ifndef SCRIPT_INTERFACE_SHAPES_SIMPLE_PORE_HPP
#define SCRIPT_INTERFACE_SHAPES_SIMPLE_PORE_HPP

#include "Shape.hpp"
#include "core/shapes/SimplePore.hpp"

namespace ScriptInterface {
namespace Shapes {

class SimplePore : public Shape {
  using CoreShape = ::Shapes::SimplePore;
  std::shared_ptr<::Shapes::SimplePore> m_simple_pore;

public:
  SimplePore() : m_simple_pore(new ::Shapes::SimplePore()) {
    add_parameters(
        {{"radius", m_simple_pore, &CoreShape::set_radius, &CoreShape::radius},
         {"length", m_simple_pore, &CoreShape::set_length, &CoreShape::length},
         {"smoothing_radius", m_simple_pore, &CoreShape::set_smoothing_radius,
          &CoreShape::smoothing_radius},
         {"axis", m_simple_pore, &CoreShape::set_axis, &CoreShape::axis},
         {"center", m_simple_pore, &CoreShape::center}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_simple_pore;
  }
};
}
}

#endif
