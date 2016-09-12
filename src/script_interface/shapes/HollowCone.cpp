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

#include "HollowCone.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap HollowCone::valid_parameters() const {
  return {{"position_x", {ParameterType::DOUBLE, true}},
          {"position_y", {ParameterType::DOUBLE, true}},
          {"position_z", {ParameterType::DOUBLE, true}},
          {"orientation_x", {ParameterType::DOUBLE, true}},
          {"orientation_y", {ParameterType::DOUBLE, true}},
          {"orientation_z", {ParameterType::DOUBLE, true}},
          {"outer_radius", {ParameterType::DOUBLE, true}},
          {"inner_radius", {ParameterType::DOUBLE, true}},
          {"width", {ParameterType::DOUBLE, true}},
          {"opening_angle", {ParameterType::DOUBLE, true}},
          {"direction", {ParameterType::DOUBLE, true}}};
}

VariantMap HollowCone::get_parameters() const {
  return {{"position_x", m_hollowcone->position_x()},
          {"position_y", m_hollowcone->position_y()},
          {"position_z", m_hollowcone->position_z()},
          {"orientation_x", m_hollowcone->orientation_x()},
          {"orientation_y", m_hollowcone->orientation_y()},
          {"orientation_z", m_hollowcone->orientation_z()},
          {"outer_radius", m_hollowcone->outer_radius()},
          {"inner_radius", m_hollowcone->inner_radius()},
          {"width", m_hollowcone->width()},
          {"opening_angle", m_hollowcone->opening_angle()},
          {"direction", m_hollowcone->direction()}};
}

void HollowCone::set_parameter(const string &name,
                               const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("position_x", m_hollowcone->position_x());
  SET_PARAMETER_HELPER("position_y", m_hollowcone->position_y());
  SET_PARAMETER_HELPER("position_z", m_hollowcone->position_z());
  SET_PARAMETER_HELPER("orientation_x", m_hollowcone->orientation_x());
  SET_PARAMETER_HELPER("orientation_y", m_hollowcone->orientation_y());
  SET_PARAMETER_HELPER("orientation_z", m_hollowcone->orientation_z());
  SET_PARAMETER_HELPER("outer_radius", m_hollowcone->outer_radius());
  SET_PARAMETER_HELPER("inner_radius", m_hollowcone->inner_radius());
  SET_PARAMETER_HELPER("width", m_hollowcone->width());
  SET_PARAMETER_HELPER("opening_angle", m_hollowcone->opening_angle());
  SET_PARAMETER_HELPER("direction", m_hollowcone->direction());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
