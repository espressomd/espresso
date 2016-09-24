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

#include "Stomatocyte.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Stomatocyte::valid_parameters() const {
  return {{"position_x", {ParameterType::DOUBLE, true}},
          {"position_y", {ParameterType::DOUBLE, true}},
          {"position_z", {ParameterType::DOUBLE, true}},
          {"orientation_x", {ParameterType::DOUBLE, true}},
          {"orientation_y", {ParameterType::DOUBLE, true}},
          {"orientation_z", {ParameterType::DOUBLE, true}},
          {"outer_radius", {ParameterType::DOUBLE, true}},
          {"inner_radius", {ParameterType::DOUBLE, true}},
          {"layer_width", {ParameterType::DOUBLE, true}},
          {"direction", {ParameterType::DOUBLE, true}}};
}

VariantMap Stomatocyte::get_parameters() const {
  return {{"position_x", m_stomatocyte->position_x()}, {"position_y", m_stomatocyte->position_y()},
          {"position_z", m_stomatocyte->position_z()}, {"orientation_x", m_stomatocyte->orientation_x()},
          {"orientation_y", m_stomatocyte->orientation_y()}, {"orientation_z", m_stomatocyte->orientation_z()}, 
          {"outer_radius", m_stomatocyte->outer_radius()}, {"inner_radius", m_stomatocyte->inner_radius()}, 
          {"layer_width", m_stomatocyte->layer_width()}, {"direction", m_stomatocyte->direction()} };
}

void Stomatocyte::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("position_x", m_stomatocyte->position_x());
  SET_PARAMETER_HELPER("position_y", m_stomatocyte->position_y());
  SET_PARAMETER_HELPER("position_z", m_stomatocyte->position_z());
  SET_PARAMETER_HELPER("orientation_x", m_stomatocyte->orientation_x());
  SET_PARAMETER_HELPER("orientation_y", m_stomatocyte->orientation_y());
  SET_PARAMETER_HELPER("orientation_z", m_stomatocyte->orientation_z());
  SET_PARAMETER_HELPER("outer_radius", m_stomatocyte->outer_radius());
  SET_PARAMETER_HELPER("inner_radius", m_stomatocyte->inner_radius());
  SET_PARAMETER_HELPER("layer_width", m_stomatocyte->layer_width());
  SET_PARAMETER_HELPER("direction", m_stomatocyte->direction());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
