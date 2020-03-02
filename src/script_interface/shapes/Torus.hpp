/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_TORUS_WALL_HPP
#define SCRIPT_INTERFACE_TORUS_WALL_HPP

#include "Shape.hpp"
#include <shapes/Torus.hpp>

namespace ScriptInterface {
namespace Shapes {

class Torus : public Shape {
  using CoreShape = ::Shapes::Torus;
  std::shared_ptr<::Shapes::Torus> m_torus;

public:
  Torus() : m_torus(new ::Shapes::Torus()) {
    add_parameters(
        {{"radius", m_torus, &CoreShape::set_radius, &CoreShape::radius},
         {"tube_radius", m_torus, &CoreShape::set_tube_radius,
          &CoreShape::tube_radius},
         {"normal", m_torus, &CoreShape::set_normal, &CoreShape::normal},
         {"center", m_torus, &CoreShape::center},
         {"direction", m_torus, &CoreShape::direction}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_torus; }
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
