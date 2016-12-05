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

#include "Rhomboid.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Rhomboid::valid_parameters() const {
  return {{"corner", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"a", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"b", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"c", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"direction", {ParameterType::DOUBLE, true}}};
}

VariantMap Rhomboid::get_parameters() const {
  return {{"corner", m_rhomboid->pos()},
          {"a", m_rhomboid->a()},
          {"b", m_rhomboid->b()},
          {"c", m_rhomboid->c()},
          {"direction", m_rhomboid->direction()}};
}

void Rhomboid::set_parameter(const string &name,
                             const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("corner", m_rhomboid->pos());
  SET_PARAMETER_HELPER("a", m_rhomboid->a());
  SET_PARAMETER_HELPER("b", m_rhomboid->b());
  SET_PARAMETER_HELPER("c", m_rhomboid->c());
  SET_PARAMETER_HELPER("direction", m_rhomboid->direction());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
