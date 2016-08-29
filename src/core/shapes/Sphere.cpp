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

#include "Sphere.hpp"
#include <cmath>

using namespace std;

#define SQR(A) ((A)*(A))

namespace Shapes {
  int Sphere::calculate_dist(const double *ppos, double *dist, double *vec)
  {
    int i;
    double fac,  c_dist;

    c_dist=0.0;
    for(i=0;i<3;i++) {
      vec[i] = pos[i] - ppos[i];
      c_dist += SQR(vec[i]);
    }
    c_dist = sqrt(c_dist);
  
    if ( direction == -1 ) {
      /* apply force towards inside the sphere */
      *dist = rad - c_dist;
      fac = *dist / c_dist;
      for(i=0;i<3;i++) vec[i] *= fac;
    } else {
      /* apply force towards outside the sphere */
      *dist = c_dist - rad;
      fac = *dist / c_dist;
      for(i=0;i<3;i++) vec[i] *= -fac;
    }
    return 0;
  }

  Parameters Sphere::all_parameters() const {
    Parameters p;
    
    p["center"] = Parameter(Parameter::Type::DOUBLE_VECTOR, 3, true);
    p["radius"] = Parameter(Parameter::Type::DOUBLE, true);
    p["direction"] = Parameter(Parameter::Type::DOUBLE, true);

    return p;
  }

std::map<std::string, Variant> Sphere::get_parameters() {
  std::map<std::string, Variant> p;

  p["center"] = pos;
  p["radius"] = rad;
  p["direction"] = direction;
    
  return p;
}

void Sphere::set_parameter(const std::string &name, const Variant &value) {
  SET_PARAMETER_HELPER("center", pos);
  SET_PARAMETER_HELPER("radius", rad);
  SET_PARAMETER_HELPER("direction", direction);
}
}

