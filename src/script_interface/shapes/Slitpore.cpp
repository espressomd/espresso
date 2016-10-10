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

#include "Slitpore.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Slitpore::valid_parameters() const {
  return {{"pore_mouth", {ParameterType::DOUBLE, true}},
          {"upper_smoothing_radius", {ParameterType::DOUBLE, true}},
          {"lower_smoothing_radius", {ParameterType::DOUBLE, true}},
          {"channel_width", {ParameterType::DOUBLE, true}}, 
          {"pore_width", {ParameterType::DOUBLE, true}},
          {"pore_length", {ParameterType::DOUBLE, true}}};
}

VariantMap Slitpore::get_parameters() const {
  return {{"pore_mouth", m_slitpore->pore_mouth()}, {"upper_smoothing_radius", m_slitpore->upper_smoothing_radius()}, 
          {"lower_smoothing_radius", m_slitpore->lower_smoothing_radius()}, 
          {"channel_width", m_slitpore->channel_width()}, 
          {"pore_width", m_slitpore->pore_width()},
          {"pore_length", m_slitpore->pore_length()}};
}

void Slitpore::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("pore_mouth", m_slitpore->pore_mouth());
  SET_PARAMETER_HELPER("upper_smoothing_radius", m_slitpore->upper_smoothing_radius());
  SET_PARAMETER_HELPER("lower_smoothing_radius", m_slitpore->lower_smoothing_radius());
  SET_PARAMETER_HELPER("channel_width", m_slitpore->channel_width());
  SET_PARAMETER_HELPER("pore_width", m_slitpore->pore_width());
  SET_PARAMETER_HELPER("pore_length", m_slitpore->pore_length());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
