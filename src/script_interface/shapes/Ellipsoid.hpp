/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_SHAPES_ELLIPSOID_HPP
#define SCRIPT_INTERFACE_SHAPES_ELLIPSOID_HPP

#include "Shape.hpp"

#include <shapes/Ellipsoid.hpp>

#include <memory>

namespace ScriptInterface {
namespace Shapes {

class Ellipsoid : public Shape {
public:
  Ellipsoid() : m_ellipsoid(std::make_shared<::Shapes::Ellipsoid>()) {
    add_parameters({{"center", m_ellipsoid->center()},
                    {"a",
                     [this](Variant const &v) {
                       m_ellipsoid->set_semiaxis_a(get_value<double>(v));
                     },
                     [this]() { return m_ellipsoid->semiaxis_a(); }},
                    {"b",
                     [this](Variant const &v) {
                       m_ellipsoid->set_semiaxis_b(get_value<double>(v));
                     },
                     [this]() { return m_ellipsoid->semiaxis_b(); }},
                    {"direction", m_ellipsoid->direction()}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_ellipsoid;
  }

private:
  std::shared_ptr<::Shapes::Ellipsoid> m_ellipsoid;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
