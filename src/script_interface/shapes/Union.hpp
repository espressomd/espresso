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

#ifndef SCRIPT_INTERFACE_SHAPES_SHAPE_UNION_HPP
#define SCRIPT_INTERFACE_SHAPES_SHAPE_UNION_HPP

#include "Shape.hpp"
#include "script_interface/ScriptObjectRegistry.hpp"
#include <shapes/Union.hpp>

namespace ScriptInterface {
namespace Shapes {

class Union : public Shape {
public:
  Union() : m_core_shape(new ::Shapes::Union()) {}

  Variant call_method(std::string const &name,
                      VariantMap const &params) override {
    if (name == "add") {
      auto const shape =
          get_value<std::shared_ptr<Shapes::Shape>>(params.at("shape"));
      m_core_shape->add(shape->shape());
      m_shapes.push_back(shape);
    } else if (name == "remove") {
      auto const shape =
          get_value<std::shared_ptr<Shapes::Shape>>(params.at("shape"));
      m_core_shape->remove(shape->shape());
      m_shapes.erase(std::remove(m_shapes.begin(), m_shapes.end(), shape),
                     m_shapes.end());
    } else if (name == "get_elements") {
      std::vector<Variant> ret;
      ret.reserve(m_shapes.size());
      for (auto const &s : m_shapes) {
        ret.emplace_back(s->id());
      }
      return ret;
    } else if (name == "clear") {
      for (auto &s : m_shapes) {
        m_core_shape->remove(s->shape());
        m_shapes.clear();
      }
    } else if (name == "size") {
      return static_cast<int>(m_shapes.size());
    } else if (name == "empty") {
      return m_shapes.empty();
    }
    return Shape::call_method(name, params);
  }
  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_core_shape;
  }

private:
  std::shared_ptr<::Shapes::Union> m_core_shape;
  std::vector<std::shared_ptr<Shapes::Shape>> m_shapes;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
