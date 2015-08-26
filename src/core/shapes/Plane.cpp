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

#include "Plane.hpp"
#include <cmath>

using namespace std;

#define SQR(A) ((A)*(A))

namespace Shapes {
  int Plane::calculate_dist(const double* ppos, double *dist, double *vec)
  {
    int i;
    double c_dist_sqr,c_dist;
  
    c_dist_sqr=0.0;
    for(i=0;i<3;i++) {
      if(pos[i] >= 0) {
        vec[i] = pos[i] - ppos[i];
        c_dist_sqr += SQR(vec[i]);
      }else{
        vec[i] = 0.0;
        c_dist_sqr += SQR(vec[i]);
      }
    }
    c_dist = sqrt(c_dist_sqr);
    *dist = c_dist;

  
    for(i=0;i<3;i++) {
      vec[i] *= -1;
    }
    return 0;
  }
}

