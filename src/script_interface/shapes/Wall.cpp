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

namespace ScriptInterface { namespace Shapes {

ParameterMap Wall::all_parameters() const {
  ParameterMap p;
  p["normal"] = Parameter(ParameterType::DOUBLE_VECTOR, 3, true);
  p["dist"] = Parameter(ParameterType::DOUBLE, true);

  return p;
}

VariantMap Wall::get_parameters() const {
  VariantMap p;

  p["normal"] = vector<double>(m_wall.n().begin(), m_wall.n().end());
  p["dist"] = m_wall.d();

  return p;
}

void Wall::set_parameter(const string &name, const ScriptInterface::Variant &value) {
  /* We need to tranform the vector<double> to a Vector3d for normal */
  if(name == "normal") {
    /* Get the variant as vector, and explicitly construct a Vector3d
       from that. */
    m_wall.set_normal(Vector3d(boost::get<vector<double> >(value).data()));   
  }

  SET_PARAMETER_HELPER("dist", m_wall.d());
}

} /* namespace Shapes */ } /* namespace ScriptInterface */

