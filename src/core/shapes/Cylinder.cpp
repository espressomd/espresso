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

#include "Cylinder.hpp"
#include <cmath>

using namespace std;

#define SQR(A) ((A)*(A))

namespace Shapes {
int Cylinder::calculate_dist(const double *ppos, double *dist, double *vec)
{
  int i;
  double d_per,d_par,d_real,d_per_vec[3],d_par_vec[3],d_real_vec[3];

  d_real = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = ppos[i] - pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);
    
  d_par=0.;
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * axis[i]);
  }
    
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * axis[i] ;
    d_per_vec[i] = ppos[i] - (pos[i] + d_par_vec[i]) ;
  }
		
  d_per=sqrt(SQR(d_real)-SQR(d_par));
  d_par = fabs(d_par) ;

  if ( direction == -1 ) {
    /*apply force towards inside cylinder */
    d_per = rad - d_per ;
    d_par = length - d_par;
    if (d_per < d_par )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
        vec[i]= -d_per_vec[i] * d_per /  (rad - d_per) ;
      }
    } else {
      *dist = d_par ;
      for (i=0; i<3;i++) {
        vec[i]= -d_par_vec[i] * d_par /  (length - d_par) ;
      }
    }
  } else {
    /*apply force towards outside cylinder */
    d_per = d_per - rad ;
    d_par = d_par - length ;
    if (d_par < 0 )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
        vec[i]= d_per_vec[i] * d_per /  (d_per + rad) ;
      }
    } else if ( d_per < 0) {
      *dist = d_par ;
      for (i=0; i<3;i++) {
        vec[i]= d_par_vec[i] * d_par /  (d_par + length) ;
      }
    } else {
      *dist = sqrt( SQR(d_par) + SQR(d_per)) ;
      for (i=0; i<3;i++) {
        vec[i]=
            d_per_vec[i] * d_per /  (d_per + rad) +
            d_par_vec[i] * d_par /  (d_par + length) ;
      }	
    }
  }
  return 0;
}

Parameters Cylinder::all_parameters() const {
  Parameters p;

  p["center"] = Parameter(Parameter::Type::DOUBLE_VECTOR, 3, true);
  p["axis"] = Parameter(Parameter::Type::DOUBLE_VECTOR, 3, true);
  p["radius"] = Parameter(Parameter::Type::DOUBLE, true);
  p["length"] = Parameter(Parameter::Type::DOUBLE, true);
  p["direction"] = Parameter(Parameter::Type::DOUBLE, true);
  
  return p;
}

std::map<std::string, Variant> Cylinder::get_parameters() {
  std::map<std::string, Variant> p;

  p["center"] = pos;
  p["axis"] = axis;
  p["radius"] = rad;
  p["length"] = length;
  p["direction"] = direction;
      
  return p;
}

void Cylinder::set_parameter(const std::string &name, const Variant &value) {
  SET_PARAMETER_HELPER("center", pos);
  SET_PARAMETER_HELPER("axis", axis);
  SET_PARAMETER_HELPER("radius", rad);
  SET_PARAMETER_HELPER("length", length);
  SET_PARAMETER_HELPER("direction", direction);
}

}
