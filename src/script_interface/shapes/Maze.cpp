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

#include "Maze.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {

ParameterMap Maze::valid_parameters() const {
  return {{"nsphere", {ParameterType::DOUBLE, true}},
          {"dim", {ParameterType::DOUBLE, true}},
          {"sphrad", {ParameterType::DOUBLE, true}},
		  {"cylrad", {ParameterType::DOUBLE, true}}};
}

VariantMap Maze::get_parameters() const {
  return {{"nsphere", m_maze->nsphere()}, {"dim", m_maze->dim()},
          {"sphrad", m_maze->sphrad()}, {"cylrad", m_maze->cylrad()}};
}

void Maze::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {

  SET_PARAMETER_HELPER("nsphere", m_maze->nsphere());
  SET_PARAMETER_HELPER("dim", m_maze->dim());
  SET_PARAMETER_HELPER("sphrad", m_maze->sphrad());
  SET_PARAMETER_HELPER("cylrad", m_maze->cylrad());
}

} /* namespace Shapes */
} /* namespace ScriptInterface */
