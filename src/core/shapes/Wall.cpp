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

#include <iostream>

using namespace std;

namespace Shapes {
Parameters Wall::all_parameters() const {
  Parameters p;
  p["normal"] = Parameter(Variant::DOUBLE_VECTOR, 3, true);
  p["dist"] = Parameter(Variant::DOUBLE, true);

  return p;
}

Parameters Wall::get_parameters() {
  Parameters p = all_parameters();

  p["normal"] = n;    
  p["dist"] = d;

  return p;
}

void Wall::set_parameter(const std::string &name, const Variant &value) {
  std::cout << "Wall::set_parameter(" << name << ", " << value << ")" << std::endl;
  SET_PARAMETER_HELPER("normal", n);
  SET_PARAMETER_HELPER("dist", d);

  /** Normalize the normal */
  n.normalize();
  
  std::cout << "n[0] = " << n[0] << std::endl;
}


int Wall::calculate_dist(const double *ppos, double *dist, double *vec)
{
  int i;

  *dist = -d;
  for(i=0;i<3;i++) *dist += ppos[i]*n[i];
  
  for(i=0;i<3;i++) vec[i] = n[i] * *dist;
  return 0;
}

}
