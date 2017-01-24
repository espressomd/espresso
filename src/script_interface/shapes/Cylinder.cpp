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

#include "Cylinder.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Cylinder::valid_parameters() const {
  return {{"center", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"axis", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"direction", {ParameterType::DOUBLE, true}},
          {"radius", {ParameterType::DOUBLE, true}},
          {"length", {ParameterType::DOUBLE, true}}};
}

VariantMap Cylinder::get_parameters() const {
  return {{"center", m_cylinder->pos()},
          {"axis", m_cylinder->axis()},
          {"direction", m_cylinder->direction()},
          {"length", m_cylinder->length()},
          {"radius", m_cylinder->rad()}};
}

void Cylinder::set_parameter(const string &name,
                             const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("center", m_cylinder->pos());
  SET_PARAMETER_HELPER("axis", m_cylinder->axis());
  SET_PARAMETER_HELPER("radius", m_cylinder->rad());
  SET_PARAMETER_HELPER("length", m_cylinder->length());
  SET_PARAMETER_HELPER("direction", m_cylinder->direction());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
