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

#include "SpheroCylinder.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap SpheroCylinder::valid_parameters() const {
  return {{"pos", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"axis", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"length", {ParameterType::DOUBLE, true}},
          {"rad", {ParameterType::DOUBLE, true}}};
}

VariantMap SpheroCylinder::get_parameters() const {
  return {{"pos", m_spherocylinder->pos()},
          {"axis", m_spherocylinder->axis()},
          {"length", m_spherocylinder->length()},
          {"rad", m_spherocylinder->rad()}};
}

void SpheroCylinder::set_parameter(const string &name,
                                   const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("pos", m_spherocylinder->pos());
  SET_PARAMETER_HELPER("axis", m_spherocylinder->axis());
  SET_PARAMETER_HELPER("length", m_spherocylinder->length());
  SET_PARAMETER_HELPER("rad", m_spherocylinder->rad());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
