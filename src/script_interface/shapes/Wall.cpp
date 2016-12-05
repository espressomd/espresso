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

#include "Wall.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Wall::valid_parameters() const {
  return {{"normal", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"dist", {ParameterType::DOUBLE, true}}};
}

VariantMap Wall::get_parameters() const {
  return {{"normal", m_wall->n()}, {"dist", m_wall->d()}};
}

void Wall::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {
  if (name == "normal") {
    /* Get the variant as vector, and explicitly construct a Vector3d
       from that. */
    m_wall->set_normal(get_value<Vector3d>(value));
  }

  SET_PARAMETER_HELPER("dist", m_wall->d());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
