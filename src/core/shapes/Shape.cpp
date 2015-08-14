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

#include "Shape.hpp"
/** @TODO: Get rid of grid.hpp for global box_l (fw) */
#include "grid.hpp"
#include <cmath>

#define SQR(A) ((A)*(A))

namespace Shapes {

  int Wall::calculate_dist(const double *ppos, double *dist, double *vec)
  {
    int i;

    *dist = -d;
    for(i=0;i<3;i++) *dist += ppos[i]*n[i];
  
    for(i=0;i<3;i++) vec[i] = n[i] * *dist;
    return 0;
  }

  int Sphere::calculate_dist(const double *ppos, double *dist, double *vec)
  {
    int i;
    double fac,  c_dist;

    c_dist=0.0;
    for(i=0;i<3;i++) {
      vec[i] = pos[i] - ppos[i];
      c_dist += SQR(vec[i]);
    }
    c_dist = std::sqrt(c_dist);
  
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

  int Rhomboid::calculate_dist(const double *ppos, double *dist, double *vec)
  {	
    double axb[3], bxc[3], axc[3];
    double A, B, C;
    double a_dot_bxc, b_dot_axc, c_dot_axb;
    double tmp;
    double d;
	
    //calculate a couple of vectors and scalars that are going to be used frequently
	
    axb[0] = a[1]*b[2] - a[2]*b[1];
    axb[1] = a[2]*b[0] - a[0]*b[2];
    axb[2] = a[0]*b[1] - a[1]*b[0];
	
    bxc[0] = b[1]*c[2] - b[2]*c[1];
    bxc[1] = b[2]*c[0] - b[0]*c[2];
    bxc[2] = b[0]*c[1] - b[1]*c[0];
	
    axc[0] = a[1]*c[2] - a[2]*c[1];
    axc[1] = a[2]*c[0] - a[0]*c[2];
    axc[2] = a[0]*c[1] - a[1]*c[0];
	
    a_dot_bxc = a[0]*bxc[0] + a[1]*bxc[1] + a[2]*bxc[2];
    b_dot_axc = b[0]*axc[0] + b[1]*axc[1] + b[2]*axc[2];
    c_dot_axb = c[0]*axb[0] + c[1]*axb[1] + c[2]*axb[2];
	
    //represent the distance from pos to ppos as a linear combination of the edge vectors.
	
    A = (ppos[0]-pos[0])*bxc[0] + (ppos[1]-pos[1])*bxc[1] + (ppos[2]-pos[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0])*axc[0] + (ppos[1]-pos[1])*axc[1] + (ppos[2]-pos[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0])*axb[0] + (ppos[1]-pos[1])*axb[1] + (ppos[2]-pos[2])*axb[2];
    C /= c_dot_axb;
	
    //the coefficients tell whether ppos lies within the cone defined by pos and the adjacent edges
	
    if(A <= 0 && B <= 0 && C <= 0)
      {
        vec[0] = ppos[0]-pos[0];
        vec[1] = ppos[1]-pos[1];
        vec[2] = ppos[2]-pos[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
        return  0;
      }
	
    //check for cone at pos+a

    A = (ppos[0]-pos[0]-a[0])*bxc[0] + (ppos[1]-pos[1]-a[1])*bxc[1] + (ppos[2]-pos[2]-a[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-a[0])*axc[0] + (ppos[1]-pos[1]-a[1])*axc[1] + (ppos[2]-pos[2]-a[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-a[0])*axb[0] + (ppos[1]-pos[1]-a[1])*axb[1] + (ppos[2]-pos[2]-a[2])*axb[2];
    C /= c_dot_axb;

    if(A >= 0 && B <= 0 && C <= 0)
      {
        vec[0] = ppos[0]-pos[0]-a[0];
        vec[1] = ppos[1]-pos[1]-a[1];
        vec[2] = ppos[2]-pos[2]-a[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
        return  0;
      }
	
    //check for cone at pos+b

    A = (ppos[0]-pos[0]-b[0])*bxc[0] + (ppos[1]-pos[1]-b[1])*bxc[1] + (ppos[2]-pos[2]-b[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-b[0])*axc[0] + (ppos[1]-pos[1]-b[1])*axc[1] + (ppos[2]-pos[2]-b[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-b[0])*axb[0] + (ppos[1]-pos[1]-b[1])*axb[1] + (ppos[2]-pos[2]-b[2])*axb[2];
    C /= c_dot_axb;

    if(A <= 0 && B >= 0 && C <= 0)
      {
        vec[0] = ppos[0]-pos[0]-b[0];
        vec[1] = ppos[1]-pos[1]-b[1];
        vec[2] = ppos[2]-pos[2]-b[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
        return  0;
      }
	
    //check for cone at pos+c

    A = (ppos[0]-pos[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-c[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-c[0])*axc[0] + (ppos[1]-pos[1]-c[1])*axc[1] + (ppos[2]-pos[2]-c[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-c[0])*axb[0] + (ppos[1]-pos[1]-c[1])*axb[1] + (ppos[2]-pos[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(A <= 0 && B <= 0 && C >= 0)
      {
        vec[0] = ppos[0]-pos[0]-c[0];
        vec[1] = ppos[1]-pos[1]-c[1];
        vec[2] = ppos[2]-pos[2]-c[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	
        return  0;
      }
	
    //check for cone at pos+a+b

    A = (ppos[0]-pos[0]-a[0]-b[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-b[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-b[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-a[0]-b[0])*axc[0] + (ppos[1]-pos[1]-a[1]-b[1])*axc[1] + (ppos[2]-pos[2]-a[2]-b[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-a[0]-b[0])*axb[0] + (ppos[1]-pos[1]-a[1]-b[1])*axb[1] + (ppos[2]-pos[2]-a[2]-b[2])*axb[2];
    C /= c_dot_axb;

    if(A >= 0 && B >= 0 && C <= 0)
      {
        vec[0] = ppos[0]-pos[0]-a[0]-b[0];
        vec[1] = ppos[1]-pos[1]-a[1]-b[1];
        vec[2] = ppos[2]-pos[2]-a[2]-b[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for cone at pos+a+c

    A = (ppos[0]-pos[0]-a[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-c[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-a[0]-c[0])*axc[0] + (ppos[1]-pos[1]-a[1]-c[1])*axc[1] + (ppos[2]-pos[2]-a[2]-c[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-a[0]-c[0])*axb[0] + (ppos[1]-pos[1]-a[1]-c[1])*axb[1] + (ppos[2]-pos[2]-a[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(A >= 0 && B <= 0 && C >= 0)
      {
        vec[0] = ppos[0]-pos[0]-a[0]-c[0];
        vec[1] = ppos[1]-pos[1]-a[1]-c[1];
        vec[2] = ppos[2]-pos[2]-a[2]-c[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for cone at pos+a+c

    A = (ppos[0]-pos[0]-b[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-b[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-b[2]-c[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-b[0]-c[0])*axc[0] + (ppos[1]-pos[1]-b[1]-c[1])*axc[1] + (ppos[2]-pos[2]-b[2]-c[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-b[0]-c[0])*axb[0] + (ppos[1]-pos[1]-b[1]-c[1])*axb[1] + (ppos[2]-pos[2]-b[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(A <= 0 && B >= 0 && C >= 0)
      {
        vec[0] = ppos[0]-pos[0]-b[0]-c[0];
        vec[1] = ppos[1]-pos[1]-b[1]-c[1];
        vec[2] = ppos[2]-pos[2]-b[2]-c[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for cone at pos+a+b+c

    A = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*bxc[2];
    A /= a_dot_bxc;	
    B = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*axc[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*axc[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*axc[2];
    B /= b_dot_axc;	
    C = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*axb[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*axb[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(A >= 0 && B >= 0 && C >= 0)
      {
        vec[0] = ppos[0]-pos[0]-a[0]-b[0]-c[0];
        vec[1] = ppos[1]-pos[1]-a[1]-b[1]-c[1];
        vec[2] = ppos[2]-pos[2]-a[2]-b[2]-c[2];
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos, a
	
    B = (ppos[0]-pos[0])*axc[0] + (ppos[1]-pos[1])*axc[1] + (ppos[2]-pos[2])*axc[2];
    B /= b_dot_axc;
    C = (ppos[0]-pos[0])*axb[0] + (ppos[1]-pos[1])*axb[1] + (ppos[2]-pos[2])*axb[2];
    C /= c_dot_axb;
	
    if(B <= 0 && C <= 0)
      {
        tmp = (ppos[0]-pos[0])*a[0] + (ppos[1]-pos[1])*a[1] + (ppos[2]-pos[2])*a[2];
        tmp /= a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
		
        vec[0] = ppos[0]-pos[0] - a[0]*tmp;
        vec[1] = ppos[1]-pos[1] - a[1]*tmp;
        vec[2] = ppos[2]-pos[2] - a[2]*tmp;
		
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
        return  0;
      }
	
    //check for prism at edge pos, b

    A = (ppos[0]-pos[0])*bxc[0] + (ppos[1]-pos[1])*bxc[1] + (ppos[2]-pos[2])*bxc[2];
    A /= a_dot_bxc;
    C = (ppos[0]-pos[0])*axb[0] + (ppos[1]-pos[1])*axb[1] + (ppos[2]-pos[2])*axb[2];
    C /= c_dot_axb;

    if(A <= 0 && C <= 0)
      {
        tmp = (ppos[0]-pos[0])*b[0] + (ppos[1]-pos[1])*b[1] + (ppos[2]-pos[2])*b[2];
        tmp /= b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
	
        vec[0] = ppos[0]-pos[0] - b[0]*tmp;
        vec[1] = ppos[1]-pos[1] - b[1]*tmp;
        vec[2] = ppos[2]-pos[2] - b[2]*tmp;
	
        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	
        return  0;
      }
	
    //check for prism at edge pos, c

    A = (ppos[0]-pos[0])*bxc[0] + (ppos[1]-pos[1])*bxc[1] + (ppos[2]-pos[2])*bxc[2];
    A /= a_dot_bxc;
    B = (ppos[0]-pos[0])*axc[0] + (ppos[1]-pos[1])*axc[1] + (ppos[2]-pos[2])*axc[2];
    B /= b_dot_axc;

    if(A <= 0 && B <= 0)
      {
        tmp = (ppos[0]-pos[0])*c[0] + (ppos[1]-pos[1])*c[1] + (ppos[2]-pos[2])*c[2];
        tmp /= c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

        vec[0] = ppos[0]-pos[0] - c[0]*tmp;
        vec[1] = ppos[1]-pos[1] - c[1]*tmp;
        vec[2] = ppos[2]-pos[2] - c[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+a, b

    A = (ppos[0]-pos[0]-a[0])*bxc[0] + (ppos[1]-pos[1]-a[1])*bxc[1] + (ppos[2]-pos[2]-a[2])*bxc[2];
    A /= a_dot_bxc;
    C = (ppos[0]-pos[0]-a[0])*axb[0] + (ppos[1]-pos[1]-a[1])*axb[1] + (ppos[2]-pos[2]-a[2])*axb[2];
    C /= c_dot_axb;

    if(A >= 0 && C <= 0)
      {
        tmp = (ppos[0]-pos[0]-a[0])*b[0] + (ppos[1]-pos[1]-a[1])*b[1] + (ppos[2]-pos[2]-a[2])*b[2];
        tmp /= b[0]*b[0] + b[1]*b[1] + b[2]*b[2];

        vec[0] = ppos[0]-pos[0]-a[0] - b[0]*tmp;
        vec[1] = ppos[1]-pos[1]-a[1] - b[1]*tmp;
        vec[2] = ppos[2]-pos[2]-a[2] - b[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+a, c

    A = (ppos[0]-pos[0]-a[0])*bxc[0] + (ppos[1]-pos[1]-a[1])*bxc[1] + (ppos[2]-pos[2]-a[2])*bxc[2];
    A /= a_dot_bxc;
    B = (ppos[0]-pos[0]-a[0])*axc[0] + (ppos[1]-pos[1]-a[1])*axc[1] + (ppos[2]-pos[2]-a[2])*axc[2];
    B /= b_dot_axc;

    if(A >= 0 && B <= 0)
      {
        tmp = (ppos[0]-pos[0]-a[0])*c[0] + (ppos[1]-pos[1]-a[1])*c[1] + (ppos[2]-pos[2]-a[2])*c[2];
        tmp /= c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

        vec[0] = ppos[0]-pos[0]-a[0] - c[0]*tmp;
        vec[1] = ppos[1]-pos[1]-a[1] - c[1]*tmp;
        vec[2] = ppos[2]-pos[2]-a[2] - c[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+b+c, c

    A = (ppos[0]-pos[0]-b[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-b[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-b[2]-c[2])*bxc[2];
    A /= a_dot_bxc;
    B = (ppos[0]-pos[0]-b[0]-c[0])*axc[0] + (ppos[1]-pos[1]-b[1]-c[1])*axc[1] + (ppos[2]-pos[2]-b[2]-c[2])*axc[2];
    B /= b_dot_axc;

    if(A <= 0 && B >= 0)
      {
        tmp = (ppos[0]-pos[0]-b[0]-c[0])*c[0] + (ppos[1]-pos[1]-b[1]-c[1])*c[1] + (ppos[2]-pos[2]-b[2]-c[2])*c[2];
        tmp /= c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

        vec[0] = ppos[0]-pos[0]-b[0]-c[0] - c[0]*tmp;
        vec[1] = ppos[1]-pos[1]-b[1]-c[1] - c[1]*tmp;
        vec[2] = ppos[2]-pos[2]-b[2]-c[2] - c[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+b+c, b

    A = (ppos[0]-pos[0]-b[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-b[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-b[2]-c[2])*bxc[2];
    A /= a_dot_bxc;
    C = (ppos[0]-pos[0]-b[0]-c[0])*axb[0] + (ppos[1]-pos[1]-b[1]-c[1])*axb[1] + (ppos[2]-pos[2]-b[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(A <= 0 && C >= 0)
      {
        tmp = (ppos[0]-pos[0]-b[0]-c[0])*b[0] + (ppos[1]-pos[1]-b[1]-c[1])*b[1] + (ppos[2]-pos[2]-b[2]-c[2])*b[2];
        tmp /= b[0]*b[0] + b[1]*b[1] + b[2]*b[2];

        vec[0] = ppos[0]-pos[0]-b[0]-c[0] - b[0]*tmp;
        vec[1] = ppos[1]-pos[1]-b[1]-c[1] - b[1]*tmp;
        vec[2] = ppos[2]-pos[2]-b[2]-c[2] - b[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+b+c, a

    B = (ppos[0]-pos[0]-b[0]-c[0])*axc[0] + (ppos[1]-pos[1]-b[1]-c[1])*axc[1] + (ppos[2]-pos[2]-b[2]-c[2])*axc[2];
    B /= b_dot_axc;
    C = (ppos[0]-pos[0]-b[0]-c[0])*axb[0] + (ppos[1]-pos[1]-b[1]-c[1])*axb[1] + (ppos[2]-pos[2]-b[2]-c[2])*axb[2];
    C /= c_dot_axb;
																
    if(B >= 0 && C >= 0)
      {
        tmp = (ppos[0]-pos[0]-b[0]-c[0])*a[0] + (ppos[1]-pos[1]-b[1]-c[1])*a[1] + (ppos[2]-pos[2]-b[2]-c[2])*a[2];
        tmp /= a[0]*a[0] + a[1]*a[1] + a[2]*a[2];

        vec[0] = ppos[0]-pos[0]-b[0]-c[0] - a[0]*tmp;
        vec[1] = ppos[1]-pos[1]-b[1]-c[1] - a[1]*tmp;
        vec[2] = ppos[2]-pos[2]-b[2]-c[2] - a[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+a+b, a

    B = (ppos[0]-pos[0]-a[0]-b[0])*axc[0] + (ppos[1]-pos[1]-a[1]-b[1])*axc[1] + (ppos[2]-pos[2]-a[2]-b[2])*axc[2];
    B /= b_dot_axc;
    C = (ppos[0]-pos[0]-a[0]-b[0])*axb[0] + (ppos[1]-pos[1]-a[1]-b[1])*axb[1] + (ppos[2]-pos[2]-a[2]-b[2])*axb[2];
    C /= c_dot_axb;

    if(B >= 0 && C <= 0)
      {
        tmp = (ppos[0]-pos[0]-a[0]-b[0])*a[0] + (ppos[1]-pos[1]-a[1]-b[1])*a[1] + (ppos[2]-pos[2]-a[2]-b[2])*a[2];
        tmp /= a[0]*a[0] + a[1]*a[1] + a[2]*a[2];

        vec[0] = ppos[0]-pos[0]-a[0]-b[0] - a[0]*tmp;
        vec[1] = ppos[1]-pos[1]-a[1]-b[1] - a[1]*tmp;
        vec[2] = ppos[2]-pos[2]-a[2]-b[2] - a[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+a+b, c

    A = (ppos[0]-pos[0]-a[0]-b[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-b[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-b[2])*bxc[2];
    A /= a_dot_bxc;
    B = (ppos[0]-pos[0]-a[0]-b[0])*axc[0] + (ppos[1]-pos[1]-a[1]-b[1])*axc[1] + (ppos[2]-pos[2]-a[2]-b[2])*axc[2];
    B /= b_dot_axc;

    if(A >= 0 && B >= 0)
      {
        tmp = (ppos[0]-pos[0]-a[0]-b[0])*c[0] + (ppos[1]-pos[1]-a[1]-b[1])*c[1] + (ppos[2]-pos[2]-a[2]-b[2])*c[2];
        tmp /= c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

        vec[0] = ppos[0]-pos[0]-a[0]-b[0] - c[0]*tmp;
        vec[1] = ppos[1]-pos[1]-a[1]-b[1] - c[1]*tmp;
        vec[2] = ppos[2]-pos[2]-a[2]-b[2] - c[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+a+c, a

    B = (ppos[0]-pos[0]-a[0]-c[0])*axc[0] + (ppos[1]-pos[1]-a[1]-c[1])*axc[1] + (ppos[2]-pos[2]-a[2]-c[2])*axc[2];
    B /= b_dot_axc;
    C = (ppos[0]-pos[0]-a[0]-c[0])*axb[0] + (ppos[1]-pos[1]-a[1]-c[1])*axb[1] + (ppos[2]-pos[2]-a[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(B <= 0 && C >= 0)
      {
        tmp = (ppos[0]-pos[0]-a[0]-c[0])*a[0] + (ppos[1]-pos[1]-a[1]-c[1])*a[1] + (ppos[2]-pos[2]-a[2]-c[2])*a[2];
        tmp /= a[0]*a[0] + a[1]*a[1] + a[2]*a[2];

        vec[0] = ppos[0]-pos[0]-a[0]-c[0] - a[0]*tmp;
        vec[1] = ppos[1]-pos[1]-a[1]-c[1] - a[1]*tmp;
        vec[2] = ppos[2]-pos[2]-a[2]-c[2] - a[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for prism at edge pos+a+c, b

    A = (ppos[0]-pos[0]-a[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-c[2])*bxc[2];
    A /= a_dot_bxc;
    C = (ppos[0]-pos[0]-a[0]-c[0])*axb[0] + (ppos[1]-pos[1]-a[1]-c[1])*axb[1] + (ppos[2]-pos[2]-a[2]-c[2])*axb[2];
    C /= c_dot_axb;

    if(A >= 0 && C >= 0)
      {
        tmp = (ppos[0]-pos[0]-a[0]-c[0])*b[0] + (ppos[1]-pos[1]-a[1]-c[1])*b[1] + (ppos[2]-pos[2]-a[2]-c[2])*b[2];
        tmp /= b[0]*b[0] + b[1]*b[1] + b[2]*b[2];

        vec[0] = ppos[0]-pos[0]-a[0]-c[0] - b[0]*tmp;
        vec[1] = ppos[1]-pos[1]-a[1]-c[1] - b[1]*tmp;
        vec[2] = ppos[2]-pos[2]-a[2]-c[2] - b[2]*tmp;

        *dist = direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        return  0;
      }
	
    //check for face with normal -axb
	
    *dist = (ppos[0]-pos[0])*axb[0] + (ppos[1]-pos[1])*axb[1] + (ppos[2]-pos[2])*axb[2];
    *dist *= -1.;
	
    if(*dist >= 0)
      {
        tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
        *dist /= tmp;
	
        vec[0] = -*dist * axb[0]/tmp;
        vec[1] = -*dist * axb[1]/tmp;
        vec[2] = -*dist * axb[2]/tmp;
	
        *dist *= direction;

        return  0;
      }
	
    //calculate distance to face with normal axc

    *dist = (ppos[0]-pos[0])*axc[0] + (ppos[1]-pos[1])*axc[1] + (ppos[2]-pos[2])*axc[2];
	
    if(*dist >= 0)
      {
        tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
        *dist /= tmp;
	
        vec[0] = *dist * axc[0]/tmp;
        vec[1] = *dist * axc[1]/tmp;
        vec[2] = *dist * axc[2]/tmp;

        *dist *= direction;

        return  0;
      }
	
    //calculate distance to face with normal -bxc

    *dist = (ppos[0]-pos[0])*bxc[0] + (ppos[1]-pos[1])*bxc[1] + (ppos[2]-pos[2])*bxc[2];
    *dist *= -1.;
	
    if(*dist >= 0)
      {
        tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
        *dist /= tmp;
		
        vec[0] = -*dist * bxc[0]/tmp;
        vec[1] = -*dist * bxc[1]/tmp;
        vec[2] = -*dist * bxc[2]/tmp;
		
        *dist *= direction;

        return  0;
      }
	
    //calculate distance to face with normal axb
	
    *dist = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*axb[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*axb[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*axb[2];
	
    if(*dist >= 0)
      {
        tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
        *dist /= tmp;

        vec[0] = *dist * axb[0]/tmp;
        vec[1] = *dist * axb[1]/tmp;
        vec[2] = *dist * axb[2]/tmp;
	
        *dist *= direction;

        return  0;
      }
																					
    //calculate distance to face with normal -axc

    *dist = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*axc[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*axc[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*axc[2];
    *dist *= -1.;
	
    if(*dist >= 0)
      {
        tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
        *dist /= tmp;

        vec[0] = -*dist * axc[0]/tmp;
        vec[1] = -*dist * axc[1]/tmp;
        vec[2] = -*dist * axc[2]/tmp;
		
        *dist *= direction;

        return  0;
      }
																		
    //calculate distance to face with normal bxc

    *dist = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*bxc[2];
	
    if(*dist >= 0)
      {
        tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
        *dist /= tmp;

        vec[0] = *dist * bxc[0]/tmp;
        vec[1] = *dist * bxc[1]/tmp;
        vec[2] = *dist * bxc[2]/tmp;
		
        *dist *= direction;

        return  0;
      }
	
    //ppos lies within rhomboid. Find nearest wall for interaction.
	 
    //check for face with normal -axb
	
    *dist = (ppos[0]-pos[0])*axb[0] + (ppos[1]-pos[1])*axb[1] + (ppos[2]-pos[2])*axb[2];
    *dist *= -1.;
    tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
    *dist /= tmp;
	
    vec[0] = -*dist * axb[0]/tmp;
    vec[1] = -*dist * axb[1]/tmp;
    vec[2] = -*dist * axb[2]/tmp;
	
    *dist *= direction;

    //calculate distance to face with normal axc

    d = (ppos[0]-pos[0])*axc[0] + (ppos[1]-pos[1])*axc[1] + (ppos[2]-pos[2])*axc[2];
    tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
    d /= tmp;
	
    if(abs(d) < abs(*dist))
      {
        vec[0] = d * axc[0]/tmp;
        vec[1] = d * axc[1]/tmp;
        vec[2] = d * axc[2]/tmp;
	
        *dist = direction * d;
      }

    //calculate distance to face with normal -bxc

    d = (ppos[0]-pos[0])*bxc[0] + (ppos[1]-pos[1])*bxc[1] + (ppos[2]-pos[2])*bxc[2];
    d *= -1.;
    tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
    d /= tmp;

    if(abs(d) < abs(*dist))
      {							
        vec[0] = -d * bxc[0]/tmp;
        vec[1] = -d * bxc[1]/tmp;
        vec[2] = -d * bxc[2]/tmp;

        *dist = direction * d;
      }
	
    //calculate distance to face with normal axb

    d = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*axb[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*axb[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*axb[2];
    tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
    d /= tmp;
	
    if(abs(d) < abs(*dist))
      {																					
        vec[0] = d * axb[0]/tmp;
        vec[1] = d * axb[1]/tmp;
        vec[2] = d * axb[2]/tmp;

        *dist = direction * d;
      }
	
    //calculate distance to face with normal -axc

    d = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*axc[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*axc[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*axc[2];
    d *= -1.;
    tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
    d /= tmp;

    if(abs(d) < abs(*dist))
      {																						
        vec[0] = -d * axc[0]/tmp;
        vec[1] = -d * axc[1]/tmp;
        vec[2] = -d * axc[2]/tmp;

        *dist = direction * d;
      }
																		
    //calculate distance to face with normal bxc

    d = (ppos[0]-pos[0]-a[0]-b[0]-c[0])*bxc[0] + (ppos[1]-pos[1]-a[1]-b[1]-c[1])*bxc[1] + (ppos[2]-pos[2]-a[2]-b[2]-c[2])*bxc[2];
    tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
    d /= tmp;

    if(abs(d) < abs(*dist))
      {																						
        vec[0] = d * bxc[0]/tmp;
        vec[1] = d * bxc[1]/tmp;
        vec[2] = d * bxc[2]/tmp;
	
        *dist = direction * d;
      }
    return 0;
  }

  static double sign(double x) {
    if (x > 0)
      return 1.;
    else
      return -1;
  }

  int Pore::calculate_dist(const double* ppos, double *dist, double *vec)
  {
    int i;
    double c_dist[3];           /* cartesian distance from pore center */
    double z , r;             /* cylindrical coordinates, coordinate system parallel to
                                 pore, origin at pore centera */
    double z_vec[3], r_vec[3]; /* cartesian vectors that correspond to these coordinates */
    double e_z[3], e_r[3];    /* unit vectors in the cylindrical coordinate system */
    /* helper variables, for performance reasons should the be move the the
     * constraint struct*/
    double slope, slope2, z_left, z_right;
    /* and another helper that is hopefully optmized out */
    double norm;
    double c1_r, c1_z, c2_r, c2_z;
    double cone_vector_r, cone_vector_z, p1_r, p1_z, dist_vector_z, dist_vector_r, temp;

     
    slope = (rad_right - rad_left)/2./(length-smoothing_radius);
    slope2 = (outer_rad_right - outer_rad_left)/2./(length-smoothing_radius);

    /* compute the position relative to the center of the pore */
    for(i=0;i<3;i++) {
      c_dist[i] = ppos[i] - pos[i];
    } 
  
    /* compute the component parallel to the pore axis */
    z =0.; 
    for(i=0;i<3;i++) {
      z += (c_dist[i] * axis[i]);
    }
  
    /* decompose the position into parallel and perpendicular to the axis */
    r = 0.;
    for(i=0;i<3;i++) {
      z_vec[i] = z * axis[i];
      r_vec[i] = c_dist[i] - z_vec[i];
      r += r_vec[i]*r_vec[i];
    }
    r = sqrt(r);


    /* calculate norm and unit vectors for both */
    norm = 0;
    for(i=0;i<3;i++) 
      norm += z_vec[i]*z_vec[i]; 
    norm = sqrt(norm);
    for(i=0;i<3;i++) 
      e_z[i] = axis[i];
    norm = 0;
    for(i=0;i<3;i++) 
      norm += r_vec[i]*r_vec[i]; 
    norm = sqrt(norm);
    for(i=0;i<3;i++) 
      e_r[i] = r_vec[i]/norm;
  
    /* c?_r/z and are the centers of the circles that are used to smooth 
     * the entrance of the pore in cylindrical coordinates*/
    c1_z = - (length - smoothing_radius);
    c2_z = + (length - smoothing_radius);
    z_left = c1_z - sign(slope) * sqrt(slope*slope/(1+slope*slope))*smoothing_radius;
    z_right = c2_z + sign(slope) * sqrt(slope*slope/(1+slope*slope))*smoothing_radius;

    c1_r = rad_left + slope * ( z_left + length ) +
      sqrt( smoothing_radius * smoothing_radius  - SQR( z_left - c1_z ) );
    c2_r = rad_left + slope * ( z_right + length ) +
      sqrt( smoothing_radius * smoothing_radius  - SQR( z_right - c2_z ) );
    c1_r = rad_left+smoothing_radius;
    c2_r = rad_right+smoothing_radius;

    double c1_or = outer_rad_left-smoothing_radius;
    double c2_or = outer_rad_right-smoothing_radius;
 
    /* Check if we are in the region of the left wall */
    if (( (r >= c1_r) && (r <= c1_or) && (z <= c1_z) )) {
      dist_vector_z=-z - length;
      dist_vector_r=0;
      *dist = -z - length;
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return 0;
    }
    /* Check if we are in the region of the right wall */
    if (( (r >= c2_r) && (r<c2_or) && (z >= c2_z) ) ) {
      dist_vector_z=-z + length;
      dist_vector_r=0;
      *dist = +z - length;
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return 0;
    }

    /* check if the particle should feel the smoothed ends or the middle of the pore */
    /* calculate aufpunkt in z direction first.   */

    /* the distance of the particle from the pore cylinder/cone calculated by projection on the
     * cone normal. Should be > 0 if particle is inside the pore */

    cone_vector_z=1/sqrt(1+slope*slope);
    cone_vector_r=slope/sqrt(1+slope*slope);
  
    double cone_vector_z_o=1/sqrt(1+slope2*slope2);
    double cone_vector_r_o=slope2/sqrt(1+slope2*slope2);

    p1_r = c1_r+ ( (r-c1_r)*cone_vector_r + (z-c1_z)*cone_vector_z) * cone_vector_r;
    p1_z = c1_z+ ( (r-c1_r)*cone_vector_r + (z-c1_z)*cone_vector_z) * cone_vector_z;

    double p2_r = c1_or+ ( (r-c1_or)*cone_vector_r_o + (z-c1_z)*cone_vector_z_o) * cone_vector_r_o;
    double p2_z = c1_z+ ( (r-c1_or)*cone_vector_r_o + (z-c1_z)*cone_vector_z_o) * cone_vector_z_o;

    dist_vector_r = p1_r-r;
    dist_vector_z = p1_z-z;

    double dist_vector_r_o = p2_r-r;
    double dist_vector_z_o = p2_z-z;

    if ( p1_z>=c1_z && p1_z<=c2_z && dist_vector_r >= 0 ) {
      temp=sqrt( dist_vector_r*dist_vector_r + dist_vector_z*dist_vector_z );
      *dist=temp-smoothing_radius;
      dist_vector_r-=dist_vector_r/temp*smoothing_radius;
      dist_vector_z-=dist_vector_z/temp*smoothing_radius;
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return 0;
    }


    if ( p2_z>=c1_z && p2_z<=c2_z && dist_vector_r_o <= 0 ) {
      temp=sqrt( dist_vector_r_o*dist_vector_r_o + dist_vector_z_o*dist_vector_z_o );
      *dist=temp-smoothing_radius;
      dist_vector_r_o-=dist_vector_r_o/temp*smoothing_radius;
      dist_vector_z_o-=dist_vector_z_o/temp*smoothing_radius;
      for (i=0; i<3; i++) vec[i]=-dist_vector_r_o*e_r[i] - dist_vector_z_o*e_z[i];
      return 0;
    }


    /* Check if we are in the range of the left smoothing circle */
    if (p1_z <= c1_z && r <= c1_r ) {
      /* distance from the smoothing center */
      norm = sqrt( (z - c1_z)*(z - c1_z) + (r - c1_r)*(r - c1_r) );
      *dist = norm - smoothing_radius;
      dist_vector_r=(smoothing_radius/norm -1)*(r - c1_r);
      dist_vector_z=(smoothing_radius/norm - 1)*(z - c1_z);
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return 0;
    }
    /* upper left smoothing circle */
    if (p2_z <= c1_z && r >= c1_or ) {
      /* distance from the smoothing center */
      norm = sqrt( (z - c1_z)*(z - c1_z) + (r - c1_or)*(r - c1_or) );
      *dist = norm - smoothing_radius;
      dist_vector_r=(smoothing_radius/norm -1)*(r - c1_or);
      dist_vector_z=(smoothing_radius/norm - 1)*(z - c1_z);
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return  0;
    }
    /* Check if we are in the range of the right smoothing circle */
    if (p1_z >= c2_z && r <= c2_r ) {
      norm = sqrt( (z - c2_z)*(z - c2_z) + (r - c2_r)*(r - c2_r) );
      *dist = norm - smoothing_radius;
      dist_vector_r=(smoothing_radius/norm -1)*(r - c2_or);
      dist_vector_z=(smoothing_radius/norm - 1)*(z - c2_z);
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return  0;
    }
    /* Check if we are in the range of the upper right smoothing circle */
    if (p2_z >= c2_z && r >= c2_or ) {
      norm = sqrt( (z - c2_z)*(z - c2_z) + (r - c2_or)*(r - c2_or) );
      *dist = norm - smoothing_radius;
      dist_vector_r=(smoothing_radius/norm -1)*(r - c2_or);
      dist_vector_z=(smoothing_radius/norm - 1)*(z - c2_z);
      for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
      return  0;
    }
    *dist=-1e99;
    vec[0] = vec[1] = vec[2] = 1e99;
    return 0;
    //  exit(printf("should never be reached, z %f, r%f\n",z, r));
  }

  int Maze::calculate_dist(const double *ppos, double *dist, double *vec)
  {
    int i,min_axis,cursph[3],dim;
    double diasph,fac,c_dist,sph_dist,cyl_dist,temp_dis;
    double sph_vec[3],cyl_vec[3];

    dim=(int) this->dim;
    diasph = box_l[0]/nsphere;

    /* First determine the distance to the sphere */
    c_dist=0.0;
    for(i=0;i<3;i++) {
      cursph[i] = (int) (ppos[i]/diasph);
      sph_vec[i] = (cursph[i]+0.5) * diasph  - (ppos[i]);
      c_dist += SQR(sph_vec[i]);
    }
    c_dist = sqrt(c_dist);
    sph_dist = sphrad - c_dist;
    fac = sph_dist / c_dist;
    for(i=0;i<3;i++) cyl_vec[i] = sph_vec[i];
    for(i=0;i<3;i++) sph_vec[i] *= fac;
  
    /* Now calculate the cylinder stuff */
    /* have to check all 3 cylinders */
    min_axis=2;
    cyl_dist=sqrt(cyl_vec[0]*cyl_vec[0]+cyl_vec[1]*cyl_vec[1]);
  
    if(dim > 0 ){
      temp_dis=sqrt(cyl_vec[0]*cyl_vec[0]+cyl_vec[2]*cyl_vec[2]);
      if ( temp_dis < cyl_dist) {
        min_axis=1;
        cyl_dist=temp_dis;
      }

      if(dim > 1 ){
        temp_dis=sqrt(cyl_vec[1]*cyl_vec[1]+cyl_vec[2]*cyl_vec[2]);
        if ( temp_dis < cyl_dist) {
          min_axis=0;
          cyl_dist=temp_dis;
        }
      }
    }
    cyl_vec[min_axis]=0.;
  
    c_dist=cyl_dist;
    cyl_dist = cylrad - c_dist;
    fac = cyl_dist / c_dist;
    for(i=0;i<3;i++) cyl_vec[i] *= fac;
  
    /* Now decide between cylinder and sphere */
    if ( sph_dist > 0 ) {
      if ( sph_dist>cyl_dist ) {
        *dist = sph_dist;
        for(i=0;i<3;i++) vec[i] = sph_vec[i];
      } else {
        *dist = cyl_dist;
        for(i=0;i<3;i++) vec[i] = cyl_vec[i];  
      }
    } else {
      *dist = cyl_dist;
      for(i=0;i<3;i++) vec[i] = cyl_vec[i];  
    }
    return 0;
  }

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

  int SpheroCylinder::calculate_dist(const double *ppos, double *dist, double *vec)
  {
    int i;
    double d = 0.0;
    double ppos_local[3];

    for(i = 0; i < 3; i++) {
      ppos_local[i] = ppos[i] - pos[i];
      d += ppos_local[i] * axis[i];
    }

    if(abs(d) >= length) {
      *dist = 0.0;
    
      for(i = 0; i < 3; i++) {
        vec[i] = ppos_local[i] - length * axis[i] * sign(d);
        *dist += vec[i]*vec[i];
      }

      *dist = sqrt(*dist);
    
      if(*dist != 0.0)
        for(i = 0; i < 3; i++)
          vec[i] /= *dist;
    
      *dist -= rad;

      for(i = 0; i < 3; i++)
        vec[i] *= *dist;

      *dist *= direction;
    }
    else
      Cylinder::calculate_dist(ppos, dist, vec);
    return 0;
  }

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

  int HollowCone::calculate_dist(const double *ppos, double *dist, double *vec) {
    int number;
    double r0, r1, w, alpha, xd, yd, zd,
      mu, x_2D, y_2D, t0, t1, t2,
      time1, time2, time3, time4,
      mdst0, mdst1, mdst2, mdst3,
      mindist, normalize, x, y, z,
      distance, normal_x, normal_y, direction,
      xp, yp, zp, xpp, ypp, sin_xy, cos_xy,
      normal_x_3D, normal_y_3D, normal_z_3D,
      normal_3D_x, normal_3D_y, normal_3D_z;

    double closest_point_3D [3] = { -1.0, -1.0, -1.0 };

    // Set the dimensions of the hollow cone

    r0 = inner_radius;
    r1 = outer_radius;
    w = width;
    alpha = opening_angle;

    // Set the position and orientation of the hollow cone

    double hollow_cone_3D_position [3] = { position_x,
                                           position_y,
                                           position_z };

    double hollow_cone_3D_orientation [3] = { orientation_x,
                                              orientation_y,
                                              orientation_z };

    // Set the point for which we want to know the distance

    double point_3D[3];
  
    point_3D[0] = ppos[0];
    point_3D[1] = ppos[1];
    point_3D[2] = ppos[2];

    /***** Convert 3D coordinates to 2D planar coordinates *****/

    // Calculate the point on position + mu * orientation,
    // where the difference segment is orthogonal

    mu = (
          hollow_cone_3D_orientation[0]*point_3D[0] 
          + hollow_cone_3D_orientation[1]*point_3D[1]  
          + hollow_cone_3D_orientation[2]*point_3D[2]
          - hollow_cone_3D_position[0]*hollow_cone_3D_orientation[0] 
          - hollow_cone_3D_position[1]*hollow_cone_3D_orientation[1] 
          - hollow_cone_3D_position[2]*hollow_cone_3D_orientation[2] 
          ) / (
               hollow_cone_3D_orientation[0]*hollow_cone_3D_orientation[0] 
               + hollow_cone_3D_orientation[1]*hollow_cone_3D_orientation[1] 
               + hollow_cone_3D_orientation[2]*hollow_cone_3D_orientation[2]
               );
  
    // Then the closest point to the line is

    closest_point_3D[0] =   hollow_cone_3D_position[0]
      + mu*hollow_cone_3D_orientation[0];
    closest_point_3D[1] =   hollow_cone_3D_position[1]
      + mu*hollow_cone_3D_orientation[1];
    closest_point_3D[2] =   hollow_cone_3D_position[2]
      + mu*hollow_cone_3D_orientation[2];

    // So the shortest distance to the line is

    x_2D = sqrt(
                ( closest_point_3D[0] - point_3D[0] ) * 
                ( closest_point_3D[0] - point_3D[0] ) 
                + ( closest_point_3D[1] - point_3D[1] ) * 
                ( closest_point_3D[1] - point_3D[1] )
                + ( closest_point_3D[2] - point_3D[2] ) * 
                ( closest_point_3D[2] - point_3D[2] )
                );

  
    y_2D = mu*sqrt(
                   hollow_cone_3D_orientation[0]*hollow_cone_3D_orientation[0] 
                   + hollow_cone_3D_orientation[1]*hollow_cone_3D_orientation[1] 
                   + hollow_cone_3D_orientation[2]*hollow_cone_3D_orientation[2]
                   );

    /***** Use the obtained planar coordinates in distance function *****/

    // Calculate intermediate results which we need to determine
    // the distance

    t0 = ( y_2D*cos(alpha) + ( x_2D - r0 )*sin(alpha) )/r1;
    t1 = ( w - 2.0*( x_2D - r0 )*cos(alpha) + 2.0*y_2D*sin(alpha) )/( 2.0*w );
    t2 = ( w + 2.0*( x_2D - r0 )*cos(alpha) - 2.0*y_2D*sin(alpha) )/( 2.0*w );

    if ( t0 >= 0.0 && t0 <= 1.0 )
      {
        time1 = t0;
        time2 = t0;
      } 
    else if ( t0 > 1.0 )
      {
        time1 = 1.0;
        time2 = 1.0;
      }
    else
      {
        time1 = 0.0;
        time2 = 0.0;
      }

    if ( t1 >= 0.0 && t1 <= 1.0 )
      {
        time3 = t1;
      } 
    else if ( t1 > 1.0 )
      {
        time3 = 1.0;
      }
    else
      {
        time3 = 0.0;
      }

    if ( t2 >= 0.0 && t2 <= 1.0 )
      {
        time4 = t2;
      } 
    else if ( t2 > 1.0 )
      {
        time4 = 1.0;
      }
    else
      {
        time4 = 0.0;
      }

    mdst0 =   x_2D*x_2D + y_2D*y_2D - 2.0*x_2D*r0 + r0*r0 + r1*r1*time1*time1
      + 0.25*w*w + ( -2.0*y_2D*r1*time1 + (  x_2D - r0 )*w )*cos(alpha)
      - (  2.0*x_2D*r1*time1 - 2.0*r0*r1*time1 + y_2D*w )*sin(alpha);
    mdst0 = sqrt(mdst0);

    mdst1 =   x_2D*x_2D + y_2D*y_2D - 2.0*x_2D*r0 + r0*r0 + r1*r1*time2*time2
      + 0.25*w*w + ( -2.0*y_2D*r1*time2 + ( -x_2D + r0 )*w )*cos(alpha)
      + ( -2.0*x_2D*r1*time2 + 2.0*r0*r1*time2 + y_2D*w )*sin(alpha);
    mdst1 = sqrt(mdst1);

    mdst2 =   x_2D*x_2D + y_2D*y_2D - 2.0*x_2D*r0 + r0*r0
      + 0.25*w*w - time3*w*w + time3*time3*w*w
      + ( x_2D - r0 )*( -1.0 + 2.0*time3 )*w*cos(alpha)
      - y_2D*( -1.0 + 2.0*time3 )*w*sin(alpha);
    mdst2 = sqrt(mdst2);

    mdst3 =   x_2D*x_2D + y_2D*y_2D - 2.0*x_2D*r0 + r0*r0
      + 0.25*w*w - time4*w*w + time4*time4*w*w
      - ( x_2D - r0 )*( -1.0 + 2.0*time4 )*w*cos(alpha)
      + y_2D*( -1.0 + 2.0*time4 )*w*sin(alpha)
      + r1*r1 - 2.0*y_2D*r1*cos(alpha)
      + ( -2.0*x_2D*r1 + 2.0*r0*r1 )*sin(alpha);
    mdst3 = sqrt(mdst3);

    double distlist[4] = { mdst0, mdst1, mdst2, mdst3 };

    // Now we only need to determine which distance is minimal
    // and remember which one it is

    mindist = -1.0;

    for ( int i = 0; i < 4; i++ )
      {
        if ( mindist < 0.0 )
          {
            number = i;
            mindist = distlist[i];
          }

        if ( mindist > distlist[i] )  
          {
            number = i;
            mindist = distlist[i];
          }
      }

    // Now we know the number corresponding to the boundary
    // to which the point is closest, we know the distance,
    // but we still need the normal

    distance = -1.0;
    normal_x = -1.0;
    normal_y = -1.0;

    if ( number == 0 )
      {
        normal_x = x_2D - r0 + 0.5*w*cos(alpha) - r1*time1*sin(alpha);
        normal_y = y_2D - r1*time1*cos(alpha) - 0.5*w*sin(alpha);
        normalize = 1.0/sqrt( normal_x*normal_x + normal_y*normal_y );
        normal_x *= normalize;
        normal_y *= normalize;

        direction = -normal_x*cos(alpha) + normal_y*sin(alpha);

        if ( fabs(direction) < 1.0e-06 && ( fabs( time1 - 0.0 ) < 1.0e-06 || fabs( time1 - 1.0 ) < 1.0e-06 ) )
          {
            if( fabs( time1 - 0.0 ) < 1.0e-06 )
              {
                direction = -normal_x*sin(alpha) - normal_y*cos(alpha);
              }
            else
              {
                direction = normal_x*sin(alpha) + normal_y*cos(alpha);
              }
          }

        if ( direction > 0.0 )
          {
            distance = mindist;
          }
        else
          {
            distance = -mindist;
            normal_x *= -1.0;
            normal_y *= -1.0; 
          }
      }
    else if ( number == 1 )
      {
        normal_x = x_2D - r0 - 0.5*w*cos(alpha) - r1*time2*sin(alpha);
        normal_y = y_2D - r1*time2*cos(alpha) + 0.5*w*sin(alpha);
        normalize = 1.0/sqrt( normal_x*normal_x + normal_y*normal_y );
        normal_x *= normalize;
        normal_y *= normalize;

        direction = normal_x*cos(alpha) - normal_y*sin(alpha);

        if ( fabs(direction) < 1.0e-06 && ( fabs( time2 - 0.0 ) < 1.0e-06 || fabs( time2 - 1.0 ) < 1.0e-06 ) )
          {
            if( fabs( time2 - 0.0 ) < 1.0e-06 )
              {
                direction = -normal_x*sin(alpha) - normal_y*cos(alpha);
              }
            else
              {
                direction = normal_x*sin(alpha) + normal_y*cos(alpha);
              }
          }

        if ( direction > 0.0 )
          {
            distance = mindist;
          }
        else
          {
            distance = -mindist;
            normal_x *= -1.0;
            normal_y *= -1.0; 
          }
      }
    else if ( number == 2 )
      {
        normal_x = x_2D - r0 - 0.5*( 1.0 - 2.0*time3 )*w*cos(alpha);
        normal_y = y_2D + 0.5*( 1.0 - 2.0*time3 )*w*sin(alpha);
        normalize = 1.0/sqrt( normal_x*normal_x + normal_y*normal_y );
        normal_x *= normalize;
        normal_y *= normalize;

        direction = -normal_x*sin(alpha) - normal_y*cos(alpha);

        if ( fabs(direction) < 1.0e-06 && ( fabs( time3 - 0.0 ) < 1.0e-06 || fabs( time3 - 1.0 ) < 1.0e-06 ) )
          {
            if( fabs( time3 - 0.0 ) < 1.0e-06 )
              {
                direction = normal_x*cos(alpha) - normal_y*sin(alpha);
              }
            else
              {
                direction = -normal_x*cos(alpha) + normal_y*sin(alpha);
              }
          }

        if ( direction > 0.0 )
          {
            distance = mindist;
          }
        else
          {
            distance = -mindist;
            normal_x *= -1.0;
            normal_y *= -1.0; 
          }
      }
    else if ( number == 3 )
      {
        normal_x = x_2D - r0 + 0.5*( 1.0 - 2.0*time4 )*w*cos(alpha) - r1*sin(alpha);
        normal_y = y_2D - 0.5*( 1.0 - 2.0*time4 )*w*sin(alpha) - r1*cos(alpha);
        normalize = 1.0/sqrt( normal_x*normal_x + normal_y*normal_y );
        normal_x *= normalize;
        normal_y *= normalize;

        direction = normal_x*sin(alpha) + normal_y*cos(alpha);

        if ( fabs(direction) < 1.0e-06 && ( fabs( time4 - 0.0 ) < 1.0e-06 || fabs( time4 - 1.0 ) < 1.0e-06 ) )
          {
            if( fabs( time4 - 0.0 ) < 1.0e-06 )
              {
                direction = -normal_x*cos(alpha) + normal_y*sin(alpha);
              }
            else
              {
                direction = normal_x*cos(alpha) - normal_y*sin(alpha);
              }
          }

        if ( direction > 0.0 )
          {
            distance = mindist;
          }
        else
          {
            distance = -mindist;
            normal_x *= -1.0;
            normal_y *= -1.0; 
          }
      }

    /***** Convert 2D normal to 3D coordinates *****/

    // Now that we have the normal in 2D we need to make a final 
    // transformation to get it in 3D. The minimum distance stays
    // the same though. We first get the normalized direction vector.

    x = hollow_cone_3D_orientation[0];
    y = hollow_cone_3D_orientation[1];
    z = hollow_cone_3D_orientation[2];

    xd = x/sqrt( x*x + y*y + z*z);
    yd = y/sqrt( x*x + y*y + z*z);
    zd = z/sqrt( x*x + y*y + z*z);

    // We now establish the rotion matrix required to go
    // form {0,0,1} to {xd,yd,zd}

    double matrix [9];

    if ( xd*xd + yd*yd > 1.e-10 ) 
      {
        // Rotation matrix
  
        matrix[0] = ( yd*yd + xd*xd*zd )/( xd*xd + yd*yd );
        matrix[1] = ( xd*yd*( zd - 1.0 ) )/( xd*xd + yd*yd );
        matrix[2] = xd;

        matrix[3] = ( xd*yd*( zd - 1.0 ) )/( xd*xd + yd*yd );
        matrix[4] = ( xd*xd + yd*yd*zd )/( xd*xd + yd*yd );
        matrix[5] = yd;

        matrix[6] = -xd;
        matrix[7] = -yd;
        matrix[8] =  zd;
      }
    else
      {
        // The matrix is the identity matrix or reverses
        // or does a 180 degree rotation to take z -> -z

        if ( zd > 0 ) 
          {
            matrix[0] = 1.0;
            matrix[1] = 0.0;
            matrix[2] = 0.0;

            matrix[3] = 0.0;
            matrix[4] = 1.0;
            matrix[5] = 0.0;

            matrix[6] = 0.0;
            matrix[7] = 0.0;
            matrix[8] = 1.0;
          }
        else
          {
            matrix[0] =  1.0;
            matrix[1] =  0.0;
            matrix[2] =  0.0;

            matrix[3] =  0.0;
            matrix[4] =  1.0;
            matrix[5] =  0.0;

            matrix[6] =  0.0;
            matrix[7] =  0.0;
            matrix[8] = -1.0;
          }
      }   

    // Next we determine the 3D vector between the center
    // of the hollow cylinder and the point of interest

    xp = point_3D[0] - hollow_cone_3D_position[0];
    yp = point_3D[1] - hollow_cone_3D_position[1];
    zp = point_3D[2] - hollow_cone_3D_position[2];

    // Now we use the inverse matrix to find the 
    // position of the point with respect to the origin
    // of the z-axis oriented hollow cone located
    // in the origin

    xpp = matrix[0]*xp + matrix[3]*yp + matrix[6]*zp;
    ypp = matrix[1]*xp + matrix[4]*yp + matrix[7]*zp;

    // Now use this direction to orient the normal

    if ( xpp*xpp + ypp*ypp > 1.0e-10 )
      {
        // The point is off the rotational symmetry
        // axis of the hollow cone

        sin_xy = ypp/sqrt( xpp*xpp + ypp*ypp );
        cos_xy = xpp/sqrt( xpp*xpp + ypp*ypp );

        normal_x_3D = cos_xy*normal_x;
        normal_y_3D = sin_xy*normal_x;
        normal_z_3D = normal_y;
      }
    else
      {
        // The point is on the rotational symmetry 
        // axis of the hollow cone; a finite distance
        // away from the center the normal might have
        // an x and y component!

        normal_x_3D = 0.0;
        normal_y_3D = 0.0;
        normal_z_3D = ( normal_y > 0.0 ? 1.0 : -1.0 );    
      }
   
    // Now we need to transform the normal back to
    // the real coordinate system

    normal_3D_x =   matrix[0]*normal_x_3D 
      + matrix[1]*normal_y_3D
      + matrix[2]*normal_z_3D;
    normal_3D_y =   matrix[3]*normal_x_3D 
      + matrix[4]*normal_y_3D
      + matrix[5]*normal_z_3D;
    normal_3D_z =   matrix[6]*normal_x_3D 
      + matrix[7]*normal_y_3D
      + matrix[8]*normal_z_3D;

    // Pass the values we obtained to ESPResSo

    if ( direction == -1 ) 
      {
        // Apply force towards inside hollow cone

        *dist = -distance;

        vec[0] = -normal_3D_x;
        vec[1] = -normal_3D_y;
        vec[2] = -normal_3D_z;
      }
    else
      {
        // Apply force towards inside hollow cone

        *dist = distance;

        vec[0] = normal_3D_x;
        vec[1] = normal_3D_y;
        vec[2] = normal_3D_z;
      }

    // And we are done with the hollow cone

    vec[0] *= *dist;
    vec[1] *= *dist;
    vec[2] *= *dist;

    return 0;
  }

  int Stomatocyte::calculate_dist(const double *ppos, double *dist, double *vec) {
    // Parameters

    int io0, io1, io2, io3, io4, number;

    double x_2D, y_2D, mu,
      T0, T1, T1p, T2, T3, T3sqrt, T3p, T4sqrt, T4,
      a, b, c, d, e,
      rad0, rad1, rad2, rad3,
      pt0x, pt0y, pt1x, pt1y, pt2x, pt2y, pt3x, pt3y,
      dst0, dst1, dst2, dst3,
      mdst0, mdst1, mdst2, mdst3, mdst4,
      t0, t1, t2, t3, t4, ttota,
      distance, mindist,
      normal_x, normal_y,
      time0, time1,
      xd, x, yd, y, zd, z,
      normal_3D_x, normal_3D_y, normal_3D_z,
      xp, yp, zp, xpp, ypp,
      normal_x_3D, normal_y_3D, normal_z_3D,
      sin_xy, cos_xy;

    double closest_point_3D [3] = { -1.0, -1.0, -1.0 };

    // Set the three dimensions of the stomatocyte

    a = outer_radius;
    b = inner_radius;
    c = layer_width;
    a = a*c;
    b = b*c;

    // Set the position and orientation of the stomatocyte 

    double stomatocyte_3D_position [3] = { position_x,
                                           position_y,
                                           position_z };

    double stomatocyte_3D_orientation [3] = { orientation_x,
                                              orientation_y,
                                              orientation_z };

    // Set the point for which we want to know the distance

    double point_3D [3];
  
    point_3D[0] = ppos[0];
    point_3D[1] = ppos[1];
    point_3D[2] = ppos[2];

    /***** Convert 3D coordinates to 2D planar coordinates *****/

    // Calculate the point on position + mu * orientation,
    // where the difference segment is orthogonal

    mu = (
          stomatocyte_3D_orientation[0]*point_3D[0] 
          + stomatocyte_3D_orientation[1]*point_3D[1]  
          + stomatocyte_3D_orientation[2]*point_3D[2]
          - stomatocyte_3D_position[0]*stomatocyte_3D_orientation[0] 
          - stomatocyte_3D_position[1]*stomatocyte_3D_orientation[1] 
          - stomatocyte_3D_position[2]*stomatocyte_3D_orientation[2] 
          ) / (
               stomatocyte_3D_orientation[0]*stomatocyte_3D_orientation[0] 
               + stomatocyte_3D_orientation[1]*stomatocyte_3D_orientation[1] 
               + stomatocyte_3D_orientation[2]*stomatocyte_3D_orientation[2]
               );
  
    // Then the closest point to the line is

    closest_point_3D[0] =   stomatocyte_3D_position[0]
      + mu*stomatocyte_3D_orientation[0];
    closest_point_3D[1] =   stomatocyte_3D_position[1]
      + mu*stomatocyte_3D_orientation[1];
    closest_point_3D[2] =   stomatocyte_3D_position[2]
      + mu*stomatocyte_3D_orientation[2];

    // So the shortest distance to the line is

    x_2D = sqrt(
                ( closest_point_3D[0] - point_3D[0] ) * 
                ( closest_point_3D[0] - point_3D[0] ) 
                + ( closest_point_3D[1] - point_3D[1] ) * 
                ( closest_point_3D[1] - point_3D[1] )
                + ( closest_point_3D[2] - point_3D[2] ) * 
                ( closest_point_3D[2] - point_3D[2] )
                );

  
    y_2D = mu*sqrt(
                   stomatocyte_3D_orientation[0]*stomatocyte_3D_orientation[0] 
                   + stomatocyte_3D_orientation[1]*stomatocyte_3D_orientation[1] 
                   + stomatocyte_3D_orientation[2]*stomatocyte_3D_orientation[2]
                   );

    /***** Use the obtained planar coordinates in distance function *****/

    // Calculate intermediate results which we need to determine
    // the distance, these are basically points where the different
    // segments of the 2D cut of the stomatocyte meet. The 
    // geometry is rather complex however and details will not be given

    d = -sqrt( b*b + 4.0*b*c - 5.0*c*c ) + a - 2.0*c - b;

    e = (   6.0*c*c - 2.0*c*d 
            + a*( -3.0*c + d )
            + sqrt( 
                   - ( a - 2.0*c )*( a - 2.0*c ) * 
                   ( -2.0*a*a + 8.0*a*c + c*c + 6.0*c*d + d*d )
                    ) 
            ) / ( 
                 2.0*( a - 2.0*c ) 
                  );

    T0 = acos(
              (   a*a*d + 5.0*c*c*d + d*d*d - a*a*e - 5.0*c*c*e 
                  + 6.0*c*d*e - 3.0*d*d*e - 6.0*c*e*e + 4.0*d*e*e - 2.0*e*e*e 
                  - ( 3.0*c + e)*sqrt( 
                                      - a*a*a*a - ( 5.0*c*c + d*d + 6.0*c*e 
                                                    - 2.0*d*e + 2.0*e*e ) *
                                      ( 5.0*c*c + d*d + 6.0*c*e
                                        - 2.0*d*e + 2.0*e*e ) 
                                      + 2.0*a*a*( 13.0*c*c + d*d + 6.0*c*e
                                                  - 2.0*d*e + 2.0*e*e )
                                       )
                  ) / (
                       2.0*a*( 9.0*c*c + d*d + 6.0*c*e - 2.0*d*e + 2.0*e*e )
                       )   
              );

    T1p = acos(
               - (   39.0*c*c*c + 3.0*c*d*d + 31.0*c*c*e
                     - 6.0*c*d*e + d*d*e + 12.0*c*e*e - 2.0*d*e*e
                     + 2.0*e*e*e - a*a*( 3.0*c + e ) 
                     - d*sqrt(
                              - a*a*a*a
                              - ( 5.0*c*c + d*d + 6.0*c*e - 2.0*d*e + 2.0*e*e) *
                              ( 5.0*c*c + d*d + 6.0*c*e - 2.0*d*e + 2.0*e*e)
                              + 2.0*a*a*(   13.0*c*c + d*d + 6.0*c*e 
                                            - 2.0*d*e + 2.0*e*e )
                              ) 
                     + e*sqrt(
                              - a*a*a*a
                              - ( 5.0*c*c + d*d + 6.0*c*e - 2.0*d*e + 2.0*e*e) *
                              ( 5.0*c*c + d*d + 6.0*c*e - 2.0*d*e + 2.0*e*e)
                              + 2.0*a*a*(   13.0*c*c + d*d + 6.0*c*e
                                            - 2.0*d*e + 2.0*e*e )
                              )
                     ) / (
                          4.0*c*( 9.0*c*c + d*d + 6.0*c*e - 2.0*d*e + 2.0*e*e )
                          )
               );

    T1 = 3.0*M_PI/4.0 - T1p;

    T2 = e;

    T3sqrt = - b*b*c*c*(   a*a*a*a - 4.0*a*a*a*( b + 2.0*c + d )
                           + 4.0*b*b*d*( 4.0*c + d ) 
                           + ( 9.0*c*c + 4.0*c*d + d*d ) * 
                           ( 9.0*c*c + 4.0*c*d + d*d)
                           + 4.0*b*( 18.0*c*c*c + 17.0*c*c*d + 6.0*c*d*d + d*d*d )
                           + 2.0*a*a*(   2.0*b*b + 17.0*c*c + 12.0*c*d 
                                         + 3.0*d*d + 6.0*b*( 2.0*c + d ) 
                                         )
                           - 4.0*a*(   18.0*c*c*c + 17.0*c*c*d + 6.0*c*d*d 
                                       + d*d*d + 2.0*b*b*( 2.0*c + d ) 
                                       + b*( 17.0*c*c + 12.0*c*d + 3.0*d*d )
                                       )
                           );

    if ( T3sqrt < 0.0 ) 
      T3sqrt = 0.0;

    T3p = acos(
               - ( - a*a*a*b + 4.0*b*b*b*c + 25.0*b*b*c*c + 34.0*b*c*c*c 
                   + 2.0*b*b*b*d + 12.0*b*b*c*d + 25.0*b*c*c*d + 3.0*b*b*d*d
                   + 6.0*b*c*d*d + b*d*d*d + 3.0*a*a*b*( b + 2.0*c + d ) 
                   - a*b*( 2.0*b*b + 25.0*c*c + 12.0*c*d 
                           + 3.0*d*d + 6.0*b*( 2.0*c + d ) )
                   + 3.0*sqrt( T3sqrt )
                   ) / (
                        4.0*b*c*(   a*a + b*b + 13.0*c*c + 4.0*c*d 
                                    + d*d + 2.0*b*( 2.0*c + d ) 
                                    - 2.0*a*( b + 2.0*c + d )
                                    )
                        )
               );

    T3 = 3.0*M_PI/4.0 - T3p;

    T4sqrt = - b*b*c*c*(   a*a*a*a - 4.0*a*a*a*( b + 2.0*c + d ) 
                           + 4.0*b*b*d*( 4.0*c + d ) 
                           + ( 9.0*c*c + 4.0*c*d + d*d ) * 
                           ( 9.0*c*c + 4.0*c*d + d*d )
                           + 4.0*b*(   18.0*c*c*c + 17.0*c*c*d
                                       + 6.0*c*d*d + d*d*d )
                           + 2.0*a*a*(   2.0*b*b + 17.0*c*c 
                                         + 12.0*c*d + 3.0*d*d 
                                         + 6.0*b*( 2.0*c + d ) )
                           - 4.0*a*(   18.0*c*c*c + 17.0*c*c*d
                                       + 6.0*c*d*d + d*d*d 
                                       + 2.0*b*b*( 2.0*c + d )
                                       + b*( 17.0*c*c + 12.0*c*d + 3.0*d*d )  
                                       )
                           );

    if ( T4sqrt < 0.0 ) 
      T4sqrt = 0.0;

    T4 = acos( 
              ( - a*a*a*b + 2.0*b*b*b*b + 8.0*b*b*b*c + 17.0*b*b*c*c
                + 18.0*b*c*c*c + 4.0*b*b*b*d + 12.0*b*b*c*d + 17.0*b*c*c*d
                + 3.0*b*b*d*d + 6.0*b*c*d*d + b*d*d*d
                + 3.0*a*a*b*( b + 2.0*c + d )
                - a*b*(   4.0*b*b + 17.0*c*c + 12.0*c*d
                          + 3.0*d*d + 6.0*b*( 2.0*c + d ) 
                          )
                - 3.0*sqrt( T4sqrt )
                ) / (
                     2.0*b*b*(   a*a + b*b + 13.0*c*c + 4.0*c*d
                                 + d*d + 2.0*b*( 2.0*c + d ) 
                                 - 2.0*a*( b + 2.0*c + d )
                                 )
                     )
               );

    // Radii for the various parts of the swimmer

    rad0 = a;
    rad1 = 2.0*c;
    rad2 = 2.0*c;
    rad3 = b;

    // Center points for the circles

    pt0x = 0.0;
    pt0y = 0.0;

    pt1x = 3.0*c + e;
    pt1y = d - e;

    pt2x = 3.0*c;
    pt2y = d;

    pt3x = 0.0;
    pt3y = a - b - 2*c;

    // Distance of point of interest to center points

    dst0 = sqrt(   ( pt0x - x_2D )*( pt0x - x_2D ) 
                   + ( pt0y - y_2D )*( pt0y - y_2D ) );
    dst1 = sqrt(   ( pt1x - x_2D )*( pt1x - x_2D )
                   + ( pt1y - y_2D )*( pt1y - y_2D ) );
    dst2 = sqrt(   ( pt2x - x_2D )*( pt2x - x_2D )
                   + ( pt2y - y_2D )*( pt2y - y_2D ) );
    dst3 = sqrt(   ( pt3x - x_2D )*( pt3x - x_2D )
                   + ( pt3y - y_2D )*( pt3y - y_2D ) );

    // Now for the minimum distances, the fourth 
    // is for the line segment

    mdst0 = ( dst0 < rad0 ? rad0 - dst0 : dst0 - rad0 );
    mdst1 = ( dst1 < rad1 ? rad1 - dst1 : dst1 - rad1 );
    mdst2 = ( dst2 < rad2 ? rad2 - dst2 : dst2 - rad2 );
    mdst3 = ( dst3 < rad3 ? rad3 - dst3 : dst3 - rad3 );

    mdst4 = fabs( (-3.0 + 2.0*sqrt( 2.0 ) )*c - d + x_2D + y_2D )/sqrt( 2.0 );

    // Now determine at which time during the parametrization
    // the minimum distance to each segment is achieved

    ttota = T0 + T1 + T2 + T3 + T4;

    t0 = acos( ( y_2D - pt0y )/sqrt(   ( x_2D - pt0x )*( x_2D - pt0x )
                                       + ( y_2D - pt0y )*( y_2D - pt0y ) ) 
               );
    t1 = acos( ( x_2D - pt1x )/sqrt(   ( x_2D - pt1x )*( x_2D - pt1x )
                                       + ( y_2D - pt1y )*( y_2D - pt1y ) )
               );
    t2 = acos( ( y_2D - pt2y )/sqrt(   ( x_2D - pt2x )*( x_2D - pt2x )
                                       + ( y_2D - pt2y )*( y_2D - pt2y ) )
               );
    t3 = acos( ( y_2D - pt3y )/sqrt(   ( x_2D - pt3x )*( x_2D - pt3x )
                                       + ( y_2D - pt3y )*( y_2D - pt3y ) )
               );

    t4 = ( 3.0*c - d + 2.0*e + 2.0*T0 + 2.0*T1 - x_2D + y_2D )/( 2.0*ttota );

    // Now we can use these times to check whether or not the 
    // point where the shortest distance is found is on the 
    // segment of the circles or line that contributes to the 
    // stomatocyte contour

    time0 = (T0 + T1)/ttota;
    time1 = (T0 + T1 + T2)/ttota;

    io0 = ( 0.0 <= t0 && t0 <= T0 ? 1 : 0 );
    io1 = ( T1p <= t1 && t1 <= 3.0*M_PI/4.0 && (y_2D <= d - e) ? 1 : 0 );
    io2 = ( T3p <= t2 && t2 <= 3.0*M_PI/4.0 && x_2D <= 3.0*c ? 1 : 0 );
    io3 = ( 0.0 <= t3 && t3 <= T4 ? 1 : 0 );

    io4 = ( time0 <= t4 && t4 <= time1 ? 1 : 0 );

    int iolist [5] = { io0, io1, io2, io3, io4 };
    double distlist [5] = { mdst0, mdst1, mdst2, mdst3, mdst4 };

    // Now we only need to consider those distances for which 
    // the io# flag is set to 1

    number = -1;
    mindist = -1.0;

    for ( int i = 0; i < 5; i++ )
      {
        if ( iolist[i] == 1 )
          {
            if ( mindist < 0.0 )
              {
                number = i;
                mindist = distlist[i];
              }

            if ( mindist > distlist[i] )  
              {
                number = i;
                mindist = distlist[i];
              }
          }
      }

    // Now we know the number corresponding to the boundary
    // to which the point is closest, we know the distance,
    // but we still need the normal

    distance = -1.0;
    normal_x = -1.0;
    normal_y = -1.0;

    if ( number == 0 )
      {
        distance = ( dst0 < rad0 ? -mindist : mindist );

        normal_x = ( x_2D - pt0x )/sqrt(   ( x_2D - pt0x )*( x_2D - pt0x )
                                           + ( y_2D - pt0y )*( y_2D - pt0y ) );
        normal_y = ( y_2D - pt0y )/sqrt(   ( x_2D - pt0x )*( x_2D - pt0x )
                                           + ( y_2D - pt0y )*( y_2D - pt0y ) );
      }
    else if ( number == 1 )
      {
        distance = ( dst1 < rad1 ? -mindist : mindist );

        normal_x = ( x_2D - pt1x )/sqrt(   ( x_2D - pt1x )*( x_2D - pt1x )
                                           + ( y_2D - pt1y )*( y_2D - pt1y ) );
        normal_y = ( y_2D - pt1y )/sqrt(   ( x_2D - pt1x )*( x_2D - pt1x )
                                           + ( y_2D - pt1y )*( y_2D - pt1y ) );
      }
    else if ( number == 2 )
      {
        distance = ( dst2 < rad2 ? -mindist : mindist );

        normal_x = ( x_2D - pt2x )/sqrt(   ( x_2D - pt2x )*( x_2D - pt2x )
                                           + ( y_2D - pt2y )*( y_2D - pt2y ) );
        normal_y = ( y_2D - pt2y )/sqrt(   ( x_2D - pt2x )*( x_2D - pt2x )
                                           + ( y_2D - pt2y )*( y_2D - pt2y ) );
      }
    else if ( number == 3 )
      {
        distance = ( dst3 < rad3 ? mindist : -mindist );

        normal_x = -( x_2D - pt3x )/sqrt(   ( x_2D - pt3x )*( x_2D - pt3x )
                                            + ( y_2D - pt3y )*( y_2D - pt3y ) );
        normal_y = -( y_2D - pt3y )/sqrt(   ( x_2D - pt3x )*( x_2D - pt3x )
                                            + ( y_2D - pt3y )*( y_2D - pt3y ) );
      }
    else if ( number == 4 )
      {
        normal_x = -1.0/sqrt(2.0);
        normal_y = normal_x;

        if ( (   a - b + c - 2.0*sqrt( 2.0 )*c
                 - sqrt( ( b - c )*( b + 5.0*c ) )
                 - x_2D - y_2D ) > 0 )
          distance = mindist;
        else
          distance = -mindist;
      }

    /***** Convert 2D normal to 3D coordinates *****/

    // Now that we have the normal in 2D we need to make a final 
    // transformation to get it in 3D. The minimum distance stays
    // the same though. We first get the normalized direction vector.

    x = stomatocyte_3D_orientation[0];
    y = stomatocyte_3D_orientation[1];
    z = stomatocyte_3D_orientation[2];

    xd = x/sqrt( x*x + y*y + z*z);
    yd = y/sqrt( x*x + y*y + z*z);
    zd = z/sqrt( x*x + y*y + z*z);

    // We now establish the rotion matrix required to go
    // form {0,0,1} to {xd,yd,zd}

    double matrix [9];

    if ( xd*xd + yd*yd > 1.e-10 ) 
      {
        // Rotation matrix
  
        matrix[0] = ( yd*yd + xd*xd*zd )/( xd*xd + yd*yd );
        matrix[1] = ( xd*yd*( zd - 1.0 ) )/( xd*xd + yd*yd );
        matrix[2] = xd;

        matrix[3] = ( xd*yd*( zd - 1.0 ) )/( xd*xd + yd*yd );
        matrix[4] = ( xd*xd + yd*yd*zd )/( xd*xd + yd*yd );
        matrix[5] = yd;

        matrix[6] = -xd;
        matrix[7] = -yd;
        matrix[8] =  zd;
      }
    else
      {
        // The matrix is the identity matrix or reverses
        // or does a 180 degree rotation to take z -> -z

        if ( zd > 0 ) 
          {
            matrix[0] = 1.0;
            matrix[1] = 0.0;
            matrix[2] = 0.0;

            matrix[3] = 0.0;
            matrix[4] = 1.0;
            matrix[5] = 0.0;

            matrix[6] = 0.0;
            matrix[7] = 0.0;
            matrix[8] = 1.0;
          }
        else
          {
            matrix[0] =  1.0;
            matrix[1] =  0.0;
            matrix[2] =  0.0;

            matrix[3] =  0.0;
            matrix[4] =  1.0;
            matrix[5] =  0.0;

            matrix[6] =  0.0;
            matrix[7] =  0.0;
            matrix[8] = -1.0;
          }
      }   

    // Next we determine the 3D vector between the center
    // of the stomatocyte and the point of interest

    xp = point_3D[0] - stomatocyte_3D_position[0];
    yp = point_3D[1] - stomatocyte_3D_position[1];
    zp = point_3D[2] - stomatocyte_3D_position[2];

    // Now we use the inverse matrix to find the 
    // position of the point with respect to the origin
    // of the z-axis oriented stomatocyte located
    // in the origin

    xpp = matrix[0]*xp + matrix[3]*yp + matrix[6]*zp;
    ypp = matrix[1]*xp + matrix[4]*yp + matrix[7]*zp;

    // Now use this direction to orient the normal

    if ( xpp*xpp + ypp*ypp > 1.e-10 )
      {
        // The point is off the rotational symmetry
        // axis of the stomatocyte

        sin_xy = ypp/sqrt( xpp*xpp + ypp*ypp );
        cos_xy = xpp/sqrt( xpp*xpp + ypp*ypp );

        normal_x_3D = cos_xy*normal_x;
        normal_y_3D = sin_xy*normal_x;
        normal_z_3D = normal_y;
      }
    else
      {
        // The point is on the rotational symmetry 
        // axis of the stomatocyte; a finite distance
        // away from the center the normal might have
        // an x and y component!

        normal_x_3D = 0.0;
        normal_y_3D = 0.0;
        normal_z_3D = ( normal_y > 0.0 ? 1.0 : -1.0 );    
      }
   
    // Now we need to transform the normal back to
    // the real coordinate system

    normal_3D_x =   matrix[0]*normal_x_3D 
      + matrix[1]*normal_y_3D
      + matrix[2]*normal_z_3D;
    normal_3D_y =   matrix[3]*normal_x_3D 
      + matrix[4]*normal_y_3D
      + matrix[5]*normal_z_3D;
    normal_3D_z =   matrix[6]*normal_x_3D 
      + matrix[7]*normal_y_3D
      + matrix[8]*normal_z_3D;

    // Pass the values we obtained to ESPResSo

    if ( direction == -1 ) 
      {
        // Apply force towards inside stomatocyte

        *dist = -distance;

        vec[0] = -normal_3D_x;
        vec[1] = -normal_3D_y;
        vec[2] = -normal_3D_z;
      }
    else
      {
        // Apply force towards inside stomatocyte

        *dist = distance;

        vec[0] = normal_3D_x;
        vec[1] = normal_3D_y;
        vec[2] = normal_3D_z;
      }

    // And we are done with the stomatocyte
  
    vec[0] *= *dist;
    vec[1] *= *dist;
    vec[2] *= *dist;

    return 0;
  }

  int Slitpore::calculate_dist(const double *ppos, double *dist, double *vec) {
    // the left circles
    double box_l_x = box_l[0];
    double c11[2] = { box_l_x/2-pore_width/2-upper_smoothing_radius, pore_mouth - upper_smoothing_radius };
    double c12[2] = { box_l_x/2-pore_width/2+lower_smoothing_radius, pore_mouth - pore_length  + lower_smoothing_radius };
    // the right circles
    double c21[2] = { box_l_x/2+pore_width/2+upper_smoothing_radius, pore_mouth - upper_smoothing_radius };
    double c22[2] = { box_l_x/2+pore_width/2-lower_smoothing_radius, pore_mouth - pore_length  + lower_smoothing_radius };

    //  printf("c11 %f %f\n", c11[0], c11[1]);
    //  printf("c12 %f %f\n", c12[0], c12[1]);
    //  printf("c21 %f %f\n", c21[0], c21[1]);
    //  printf("c22 %f %f\n", c22[0], c22[1]);


    if (ppos[2] > pore_mouth + channel_width/2) {
      //    printf("upper wall\n");
      // Feel the upper wall
      *dist = pore_mouth + channel_width - ppos[2];
      vec[0] = vec[1] = 0;
      vec[2] = -*dist;
      return 0;
    }

    if (ppos[0]<c11[0] || ppos[0] > c21[0]) {
      // Feel the lower wall of the channel
      //    printf("lower wall\n");
      *dist = ppos[2] - pore_mouth;
      vec[0] = vec[1] = 0;
      vec[2] = *dist;
      return 0;
    }

    if (ppos[2] > c11[1]) {
      // Feel the upper smoothing
      if (ppos[0] < box_l_x/2) {
        //    printf("upper smoothing left\n");
        *dist = sqrt( SQR(c11[0] - ppos[0]) + SQR(c11[1] - ppos[2])) - upper_smoothing_radius;
        vec[0] = -( c11[0] - ppos[0] ) * (*dist)/(*dist+upper_smoothing_radius);
        vec[1] = 0;
        vec[2] = -( c11[1] - ppos[2] ) * (*dist)/(*dist+upper_smoothing_radius);
        return 0;
      } else {
        //    printf("upper smoothing right\n");
        *dist = sqrt( SQR(c21[0] - ppos[0]) + SQR(c21[1] - ppos[2])) - upper_smoothing_radius;
        vec[0] = -( c21[0] - ppos[0] ) * (*dist)/(*dist+upper_smoothing_radius);
        vec[1] = 0;
        vec[2] = -( c21[1] - ppos[2] ) * (*dist)/(*dist+upper_smoothing_radius);
        return 0;
      }
    }
  
    if (ppos[2] > c12[1]) {
      // Feel the pore wall
      if (ppos[0] < box_l_x/2) {
        //    printf("pore left\n");
        *dist = ppos[0] - (box_l_x/2-pore_width/2);
        vec[0]=*dist;
        vec[1]=vec[2]=0;
        return 0;
      } else {
        //    printf("pore right\n");
        *dist =  (box_l_x/2+pore_width/2) - ppos[0];
        vec[0]=-*dist;
        vec[1]=vec[2]=0;
        return 0;
      }
    }

    if (ppos[0]>c12[0] && ppos[0] < c22[0]) {
      //    printf("pore end\n");
      // Feel the pore end wall
      *dist = ppos[2] - (pore_mouth-pore_length);
      vec[0]=vec[1]=0;
      vec[2]=*dist;
      return 0;
    }
    // Else
    // Feel the lower smoothing
    if (ppos[0] < box_l_x/2) {
      //    printf("lower smoothing left\n");
      *dist = -sqrt( SQR(c12[0] - ppos[0]) + SQR(c12[1] - ppos[2])) + lower_smoothing_radius;
      vec[0] = ( c12[0] - ppos[0] ) * (*dist)/(-*dist+lower_smoothing_radius);
      vec[1] = 0;
      vec[2] = ( c12[1] - ppos[2] ) * (*dist)/(-*dist+lower_smoothing_radius);
      return 0;
    } else {
      //    printf("lower smoothing right\n");
      *dist = -sqrt( SQR(c22[0] - ppos[0]) + SQR(c22[1] - ppos[2])) + lower_smoothing_radius;
      vec[0] = ( c22[0] - ppos[0] ) * (*dist)/(-*dist+lower_smoothing_radius);
      vec[1] = 0;
      vec[2] = ( c22[1] - ppos[2] ) * (*dist)/(-*dist+lower_smoothing_radius);
      return 0;
    }

    return 0;
  }
}
