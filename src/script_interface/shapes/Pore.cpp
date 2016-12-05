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

#include "Pore.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Pore::valid_parameters() const {
  return {{"pos", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"axis", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"rad_left", {ParameterType::DOUBLE, true}},
          {"rad_right", {ParameterType::DOUBLE, true}},
          {"smoothing_radius", {ParameterType::DOUBLE, true}},
          {"length", {ParameterType::DOUBLE, true}},
          {"outer_rad_left", {ParameterType::DOUBLE, true}},
          {"outer_rad_right", {ParameterType::DOUBLE, true}}};
}

VariantMap Pore::get_parameters() const {
  return {{"pos", m_pore->pos()},
          {"axis", m_pore->axis()},
          {"rad_left", m_pore->rad_left()},
          {"rad_right", m_pore->rad_right()},
          {"smoothing_radius", m_pore->smoothing_radius()},
          {"length", m_pore->smoothing_radius()},
          {"outer_rad_left", m_pore->outer_rad_left()},
          {"outer_rad_right", m_pore->outer_rad_right()}};
}

void Pore::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("pos", m_pore->pos());
  SET_PARAMETER_HELPER("axis", m_pore->axis());
  SET_PARAMETER_HELPER("rad_left", m_pore->rad_left());
  SET_PARAMETER_HELPER("rad_right", m_pore->rad_right());
  SET_PARAMETER_HELPER("smoothing_radius", m_pore->smoothing_radius());
  SET_PARAMETER_HELPER("length", m_pore->length());
  SET_PARAMETER_HELPER("outer_rad_left", m_pore->outer_rad_left());
  SET_PARAMETER_HELPER("outer_rad_right", m_pore->outer_rad_right());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
