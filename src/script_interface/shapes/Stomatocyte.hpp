/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef SCRIPT_INTERFACE_SHAPES_STOMATOCYTE_HPP
#define SCRIPT_INTERFACE_SHAPES_STOMATOCYTE_HPP

#include "Shape.hpp"
#include "core/shapes/Stomatocyte.hpp"

namespace ScriptInterface {
namespace Shapes {

class Stomatocyte : public Shape {
public:
  Stomatocyte() : m_stomatocyte(new ::Shapes::Stomatocyte()) {
    add_parameters({{"position_x", m_stomatocyte->position_x()},
         {"position_y", m_stomatocyte->position_y()},
         {"position_z", m_stomatocyte->position_z()},
         {"orientation_x", m_stomatocyte->orientation_x()},
         {"orientation_y", m_stomatocyte->orientation_y()},
         {"orientation_z", m_stomatocyte->orientation_z()},
         {"outer_radius", m_stomatocyte->outer_radius()},
         {"inner_radius", m_stomatocyte->inner_radius()},
         {"layer_width", m_stomatocyte->layer_width()},
         {"direction", m_stomatocyte->direction()}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_stomatocyte;
  }

private:
  std::shared_ptr<::Shapes::Stomatocyte> m_stomatocyte;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
