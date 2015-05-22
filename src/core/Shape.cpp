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
#include <cmath>

#define SQR(A) ((A)*(A))

int Wall::calculate_dist(double *ppos, double *dist, double *vec)
{
  int i;

  *dist = -this->d;
  for(i=0;i<3;i++) *dist += ppos[i]*this->n[i];
  
  for(i=0;i<3;i++) vec[i] = this->n[i] * *dist;
  return 0;
}

int Sphere::calculate_dist(double *ppos, double *dist, double *vec)
{
  int i;
  double fac,  c_dist;

  c_dist=0.0;
  for(i=0;i<3;i++) {
    vec[i] = this->pos[i] - ppos[i];
    c_dist += SQR(vec[i]);
  }
  c_dist = std::sqrt(c_dist);
  
  if ( this->direction == -1 ) {
    /* apply force towards inside the sphere */
    *dist = this->rad - c_dist;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= fac;
  } else {
    /* apply force towards outside the sphere */
    *dist = c_dist - this->rad;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= -fac;
  }
  return 0;
}

int Rhomboid::calculate_dist(double *ppos, double *dist, double *vec)
{	
  double axb[3], bxc[3], axc[3];
  double A, B, C;
  double a_dot_bxc, b_dot_axc, c_dot_axb;
  double tmp;
  double d;
	
  //calculate a couple of vectors and scalars that are going to be used frequently
	
  axb[0] = this->a[1]*this->b[2] - this->a[2]*this->b[1];
  axb[1] = this->a[2]*this->b[0] - this->a[0]*this->b[2];
  axb[2] = this->a[0]*this->b[1] - this->a[1]*this->b[0];
	
  bxc[0] = this->b[1]*this->c[2] - this->b[2]*this->c[1];
  bxc[1] = this->b[2]*this->c[0] - this->b[0]*this->c[2];
  bxc[2] = this->b[0]*this->c[1] - this->b[1]*this->c[0];
	
  axc[0] = this->a[1]*this->c[2] - this->a[2]*this->c[1];
  axc[1] = this->a[2]*this->c[0] - this->a[0]*this->c[2];
  axc[2] = this->a[0]*this->c[1] - this->a[1]*this->c[0];
	
  a_dot_bxc = this->a[0]*bxc[0] + this->a[1]*bxc[1] + this->a[2]*bxc[2];
  b_dot_axc = this->b[0]*axc[0] + this->b[1]*axc[1] + this->b[2]*axc[2];
  c_dot_axb = this->c[0]*axb[0] + this->c[1]*axb[1] + this->c[2]*axb[2];
	
  //represent the distance from pos to ppos as a linear combination of the edge vectors.
	
  A = (ppos[0]-this->pos[0])*bxc[0] + (ppos[1]-this->pos[1])*bxc[1] + (ppos[2]-this->pos[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0])*axc[0] + (ppos[1]-this->pos[1])*axc[1] + (ppos[2]-this->pos[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0])*axb[0] + (ppos[1]-this->pos[1])*axb[1] + (ppos[2]-this->pos[2])*axb[2];
  C /= c_dot_axb;
	
  //the coefficients tell whether ppos lies within the cone defined by pos and the adjacent edges
	
  if(A <= 0 && B <= 0 && C <= 0)
    {
      vec[0] = ppos[0]-this->pos[0];
      vec[1] = ppos[1]-this->pos[1];
      vec[2] = ppos[2]-this->pos[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
      return  0;
    }
	
  //check for cone at pos+a

  A = (ppos[0]-this->pos[0]-this->a[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->a[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->a[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2])*axb[2];
  C /= c_dot_axb;

  if(A >= 0 && B <= 0 && C <= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->a[0];
      vec[1] = ppos[1]-this->pos[1]-this->a[1];
      vec[2] = ppos[2]-this->pos[2]-this->a[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
      return  0;
    }
	
  //check for cone at pos+b

  A = (ppos[0]-this->pos[0]-this->b[0])*bxc[0] + (ppos[1]-this->pos[1]-this->b[1])*bxc[1] + (ppos[2]-this->pos[2]-this->b[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->b[0])*axc[0] + (ppos[1]-this->pos[1]-this->b[1])*axc[1] + (ppos[2]-this->pos[2]-this->b[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->b[0])*axb[0] + (ppos[1]-this->pos[1]-this->b[1])*axb[1] + (ppos[2]-this->pos[2]-this->b[2])*axb[2];
  C /= c_dot_axb;

  if(A <= 0 && B >= 0 && C <= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->b[0];
      vec[1] = ppos[1]-this->pos[1]-this->b[1];
      vec[2] = ppos[2]-this->pos[2]-this->b[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
      return  0;
    }
	
  //check for cone at pos+c

  A = (ppos[0]-this->pos[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->c[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(A <= 0 && B <= 0 && C >= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->c[0];
      vec[1] = ppos[1]-this->pos[1]-this->c[1];
      vec[2] = ppos[2]-this->pos[2]-this->c[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	
      return  0;
    }
	
  //check for cone at pos+a+b

  A = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*axb[2];
  C /= c_dot_axb;

  if(A >= 0 && B >= 0 && C <= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->b[0];
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->b[1];
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->b[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for cone at pos+a+c

  A = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(A >= 0 && B <= 0 && C >= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->c[0];
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->c[1];
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->c[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for cone at pos+a+c

  A = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(A <= 0 && B >= 0 && C >= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->b[0]-this->c[0];
      vec[1] = ppos[1]-this->pos[1]-this->b[1]-this->c[1];
      vec[2] = ppos[2]-this->pos[2]-this->b[2]-this->c[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for cone at pos+a+b+c

  A = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;	
  B = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*axc[2];
  B /= b_dot_axc;	
  C = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(A >= 0 && B >= 0 && C >= 0)
    {
      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0];
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1];
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2];
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos, a
	
  B = (ppos[0]-this->pos[0])*axc[0] + (ppos[1]-this->pos[1])*axc[1] + (ppos[2]-this->pos[2])*axc[2];
  B /= b_dot_axc;
  C = (ppos[0]-this->pos[0])*axb[0] + (ppos[1]-this->pos[1])*axb[1] + (ppos[2]-this->pos[2])*axb[2];
  C /= c_dot_axb;
	
  if(B <= 0 && C <= 0)
    {
      tmp = (ppos[0]-this->pos[0])*this->a[0] + (ppos[1]-this->pos[1])*this->a[1] + (ppos[2]-this->pos[2])*this->a[2];
      tmp /= this->a[0]*this->a[0] + this->a[1]*this->a[1] + this->a[2]*this->a[2];
		
      vec[0] = ppos[0]-this->pos[0] - this->a[0]*tmp;
      vec[1] = ppos[1]-this->pos[1] - this->a[1]*tmp;
      vec[2] = ppos[2]-this->pos[2] - this->a[2]*tmp;
		
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
      return  0;
    }
	
  //check for prism at edge pos, b

  A = (ppos[0]-this->pos[0])*bxc[0] + (ppos[1]-this->pos[1])*bxc[1] + (ppos[2]-this->pos[2])*bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0]-this->pos[0])*axb[0] + (ppos[1]-this->pos[1])*axb[1] + (ppos[2]-this->pos[2])*axb[2];
  C /= c_dot_axb;

  if(A <= 0 && C <= 0)
    {
      tmp = (ppos[0]-this->pos[0])*this->b[0] + (ppos[1]-this->pos[1])*this->b[1] + (ppos[2]-this->pos[2])*this->b[2];
      tmp /= this->b[0]*this->b[0] + this->b[1]*this->b[1] + this->b[2]*this->b[2];
	
      vec[0] = ppos[0]-this->pos[0] - this->b[0]*tmp;
      vec[1] = ppos[1]-this->pos[1] - this->b[1]*tmp;
      vec[2] = ppos[2]-this->pos[2] - this->b[2]*tmp;
	
      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	
      return  0;
    }
	
  //check for prism at edge pos, c

  A = (ppos[0]-this->pos[0])*bxc[0] + (ppos[1]-this->pos[1])*bxc[1] + (ppos[2]-this->pos[2])*bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0]-this->pos[0])*axc[0] + (ppos[1]-this->pos[1])*axc[1] + (ppos[2]-this->pos[2])*axc[2];
  B /= b_dot_axc;

  if(A <= 0 && B <= 0)
    {
      tmp = (ppos[0]-this->pos[0])*this->c[0] + (ppos[1]-this->pos[1])*this->c[1] + (ppos[2]-this->pos[2])*this->c[2];
      tmp /= this->c[0]*this->c[0] + this->c[1]*this->c[1] + this->c[2]*this->c[2];

      vec[0] = ppos[0]-this->pos[0] - this->c[0]*tmp;
      vec[1] = ppos[1]-this->pos[1] - this->c[1]*tmp;
      vec[2] = ppos[2]-this->pos[2] - this->c[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+a, b

  A = (ppos[0]-this->pos[0]-this->a[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2])*bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0]-this->pos[0]-this->a[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2])*axb[2];
  C /= c_dot_axb;

  if(A >= 0 && C <= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->a[0])*this->b[0] + (ppos[1]-this->pos[1]-this->a[1])*this->b[1] + (ppos[2]-this->pos[2]-this->a[2])*this->b[2];
      tmp /= this->b[0]*this->b[0] + this->b[1]*this->b[1] + this->b[2]*this->b[2];

      vec[0] = ppos[0]-this->pos[0]-this->a[0] - this->b[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->a[1] - this->b[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->a[2] - this->b[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+a, c

  A = (ppos[0]-this->pos[0]-this->a[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2])*bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0]-this->pos[0]-this->a[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2])*axc[2];
  B /= b_dot_axc;

  if(A >= 0 && B <= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->a[0])*this->c[0] + (ppos[1]-this->pos[1]-this->a[1])*this->c[1] + (ppos[2]-this->pos[2]-this->a[2])*this->c[2];
      tmp /= this->c[0]*this->c[0] + this->c[1]*this->c[1] + this->c[2]*this->c[2];

      vec[0] = ppos[0]-this->pos[0]-this->a[0] - this->c[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->a[1] - this->c[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->a[2] - this->c[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+b+c, c

  A = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*axc[2];
  B /= b_dot_axc;

  if(A <= 0 && B >= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*this->c[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*this->c[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*this->c[2];
      tmp /= this->c[0]*this->c[0] + this->c[1]*this->c[1] + this->c[2]*this->c[2];

      vec[0] = ppos[0]-this->pos[0]-this->b[0]-this->c[0] - this->c[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->b[1]-this->c[1] - this->c[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->b[2]-this->c[2] - this->c[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+b+c, b

  A = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(A <= 0 && C >= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*this->b[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*this->b[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*this->b[2];
      tmp /= this->b[0]*this->b[0] + this->b[1]*this->b[1] + this->b[2]*this->b[2];

      vec[0] = ppos[0]-this->pos[0]-this->b[0]-this->c[0] - this->b[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->b[1]-this->c[1] - this->b[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->b[2]-this->c[2] - this->b[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+b+c, a

  B = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*axc[2];
  B /= b_dot_axc;
  C = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*axb[2];
  C /= c_dot_axb;
																
  if(B >= 0 && C >= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->b[0]-this->c[0])*this->a[0] + (ppos[1]-this->pos[1]-this->b[1]-this->c[1])*this->a[1] + (ppos[2]-this->pos[2]-this->b[2]-this->c[2])*this->a[2];
      tmp /= this->a[0]*this->a[0] + this->a[1]*this->a[1] + this->a[2]*this->a[2];

      vec[0] = ppos[0]-this->pos[0]-this->b[0]-this->c[0] - this->a[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->b[1]-this->c[1] - this->a[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->b[2]-this->c[2] - this->a[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+a+b, a

  B = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*axc[2];
  B /= b_dot_axc;
  C = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*axb[2];
  C /= c_dot_axb;

  if(B >= 0 && C <= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*this->a[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*this->a[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*this->a[2];
      tmp /= this->a[0]*this->a[0] + this->a[1]*this->a[1] + this->a[2]*this->a[2];

      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->b[0] - this->a[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->b[1] - this->a[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->b[2] - this->a[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+a+b, c

  A = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*axc[2];
  B /= b_dot_axc;

  if(A >= 0 && B >= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->a[0]-this->b[0])*this->c[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1])*this->c[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2])*this->c[2];
      tmp /= this->c[0]*this->c[0] + this->c[1]*this->c[1] + this->c[2]*this->c[2];

      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->b[0] - this->c[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->b[1] - this->c[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->b[2] - this->c[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+a+c, a

  B = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*axc[2];
  B /= b_dot_axc;
  C = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(B <= 0 && C >= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*this->a[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*this->a[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*this->a[2];
      tmp /= this->a[0]*this->a[0] + this->a[1]*this->a[1] + this->a[2]*this->a[2];

      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->c[0] - this->a[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->c[1] - this->a[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->c[2] - this->a[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for prism at edge pos+a+c, b

  A = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*axb[2];
  C /= c_dot_axb;

  if(A >= 0 && C >= 0)
    {
      tmp = (ppos[0]-this->pos[0]-this->a[0]-this->c[0])*this->b[0] + (ppos[1]-this->pos[1]-this->a[1]-this->c[1])*this->b[1] + (ppos[2]-this->pos[2]-this->a[2]-this->c[2])*this->b[2];
      tmp /= this->b[0]*this->b[0] + this->b[1]*this->b[1] + this->b[2]*this->b[2];

      vec[0] = ppos[0]-this->pos[0]-this->a[0]-this->c[0] - this->b[0]*tmp;
      vec[1] = ppos[1]-this->pos[1]-this->a[1]-this->c[1] - this->b[1]*tmp;
      vec[2] = ppos[2]-this->pos[2]-this->a[2]-this->c[2] - this->b[2]*tmp;

      *dist = this->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      return  0;
    }
	
  //check for face with normal -axb
	
  *dist = (ppos[0]-this->pos[0])*axb[0] + (ppos[1]-this->pos[1])*axb[1] + (ppos[2]-this->pos[2])*axb[2];
  *dist *= -1.;
	
  if(*dist >= 0)
    {
      tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
      *dist /= tmp;
	
      vec[0] = -*dist * axb[0]/tmp;
      vec[1] = -*dist * axb[1]/tmp;
      vec[2] = -*dist * axb[2]/tmp;
	
      *dist *= this->direction;

      return  0;
    }
	
  //calculate distance to face with normal axc

  *dist = (ppos[0]-this->pos[0])*axc[0] + (ppos[1]-this->pos[1])*axc[1] + (ppos[2]-this->pos[2])*axc[2];
	
  if(*dist >= 0)
    {
      tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
      *dist /= tmp;
	
      vec[0] = *dist * axc[0]/tmp;
      vec[1] = *dist * axc[1]/tmp;
      vec[2] = *dist * axc[2]/tmp;

      *dist *= this->direction;

      return  0;
    }
	
  //calculate distance to face with normal -bxc

  *dist = (ppos[0]-this->pos[0])*bxc[0] + (ppos[1]-this->pos[1])*bxc[1] + (ppos[2]-this->pos[2])*bxc[2];
  *dist *= -1.;
	
  if(*dist >= 0)
    {
      tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
      *dist /= tmp;
		
      vec[0] = -*dist * bxc[0]/tmp;
      vec[1] = -*dist * bxc[1]/tmp;
      vec[2] = -*dist * bxc[2]/tmp;
		
      *dist *= this->direction;

      return  0;
    }
	
  //calculate distance to face with normal axb
	
  *dist = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*axb[2];
	
  if(*dist >= 0)
    {
      tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
      *dist /= tmp;

      vec[0] = *dist * axb[0]/tmp;
      vec[1] = *dist * axb[1]/tmp;
      vec[2] = *dist * axb[2]/tmp;
	
      *dist *= this->direction;

      return  0;
    }
																					
  //calculate distance to face with normal -axc

  *dist = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*axc[2];
  *dist *= -1.;
	
  if(*dist >= 0)
    {
      tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
      *dist /= tmp;

      vec[0] = -*dist * axc[0]/tmp;
      vec[1] = -*dist * axc[1]/tmp;
      vec[2] = -*dist * axc[2]/tmp;
		
      *dist *= this->direction;

      return  0;
    }
																		
  //calculate distance to face with normal bxc

  *dist = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*bxc[2];
	
  if(*dist >= 0)
    {
      tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
      *dist /= tmp;

      vec[0] = *dist * bxc[0]/tmp;
      vec[1] = *dist * bxc[1]/tmp;
      vec[2] = *dist * bxc[2]/tmp;
		
      *dist *= this->direction;

      return  0;
    }
	
  //ppos lies within rhomboid. Find nearest wall for interaction.
	 
  //check for face with normal -axb
	
  *dist = (ppos[0]-this->pos[0])*axb[0] + (ppos[1]-this->pos[1])*axb[1] + (ppos[2]-this->pos[2])*axb[2];
  *dist *= -1.;
  tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
  *dist /= tmp;
	
  vec[0] = -*dist * axb[0]/tmp;
  vec[1] = -*dist * axb[1]/tmp;
  vec[2] = -*dist * axb[2]/tmp;
	
  *dist *= this->direction;

  //calculate distance to face with normal axc

  d = (ppos[0]-this->pos[0])*axc[0] + (ppos[1]-this->pos[1])*axc[1] + (ppos[2]-this->pos[2])*axc[2];
  tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
  d /= tmp;
	
  if(abs(d) < abs(*dist))
    {
      vec[0] = d * axc[0]/tmp;
      vec[1] = d * axc[1]/tmp;
      vec[2] = d * axc[2]/tmp;
	
      *dist = this->direction * d;
    }

  //calculate distance to face with normal -bxc

  d = (ppos[0]-this->pos[0])*bxc[0] + (ppos[1]-this->pos[1])*bxc[1] + (ppos[2]-this->pos[2])*bxc[2];
  d *= -1.;
  tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
  d /= tmp;

  if(abs(d) < abs(*dist))
    {							
      vec[0] = -d * bxc[0]/tmp;
      vec[1] = -d * bxc[1]/tmp;
      vec[2] = -d * bxc[2]/tmp;

      *dist = this->direction * d;
    }
	
  //calculate distance to face with normal axb

  d = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*axb[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*axb[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*axb[2];
  tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
  d /= tmp;
	
  if(abs(d) < abs(*dist))
    {																					
      vec[0] = d * axb[0]/tmp;
      vec[1] = d * axb[1]/tmp;
      vec[2] = d * axb[2]/tmp;

      *dist = this->direction * d;
    }
	
  //calculate distance to face with normal -axc

  d = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*axc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*axc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*axc[2];
  d *= -1.;
  tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
  d /= tmp;

  if(abs(d) < abs(*dist))
    {																						
      vec[0] = -d * axc[0]/tmp;
      vec[1] = -d * axc[1]/tmp;
      vec[2] = -d * axc[2]/tmp;

      *dist = this->direction * d;
    }
																		
  //calculate distance to face with normal bxc

  d = (ppos[0]-this->pos[0]-this->a[0]-this->b[0]-this->c[0])*bxc[0] + (ppos[1]-this->pos[1]-this->a[1]-this->b[1]-this->c[1])*bxc[1] + (ppos[2]-this->pos[2]-this->a[2]-this->b[2]-this->c[2])*bxc[2];
  tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
  d /= tmp;

  if(abs(d) < abs(*dist))
    {																						
      vec[0] = d * bxc[0]/tmp;
      vec[1] = d * bxc[1]/tmp;
      vec[2] = d * bxc[2]/tmp;
	
      *dist = this->direction * d;
    }
  return 0;
}

static double sign(double x) {
  if (x > 0)
    return 1.;
  else
    return -1;
}

int Pore::calculate_dist(double* ppos, double *dist, double *vec)
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

     
     slope = (this->rad_right - this->rad_left)/2./(this->length-this->smoothing_radius);
     slope2 = (this->outer_rad_right - this->outer_rad_left)/2./(this->length-this->smoothing_radius);

  /* compute the position relative to the center of the pore */
  for(i=0;i<3;i++) {
    c_dist[i] = ppos[i] - this->pos[i];
  } 
  
  /* compute the component parallel to the pore axis */
  z =0.; 
  for(i=0;i<3;i++) {
    z += (c_dist[i] * this->axis[i]);
  }
  
  /* decompose the position into parallel and perpendicular to the axis */
  r = 0.;
  for(i=0;i<3;i++) {
    z_vec[i] = z * this->axis[i];
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
    e_z[i] = this->axis[i];
  norm = 0;
  for(i=0;i<3;i++) 
    norm += r_vec[i]*r_vec[i]; 
  norm = sqrt(norm);
  for(i=0;i<3;i++) 
    e_r[i] = r_vec[i]/norm;
  
  /* c?_r/z and are the centers of the circles that are used to smooth 
   * the entrance of the pore in cylindrical coordinates*/
  c1_z = - (this->length - this->smoothing_radius);
  c2_z = + (this->length - this->smoothing_radius);
  z_left = c1_z - sign(slope) * sqrt(slope*slope/(1+slope*slope))*this->smoothing_radius;
  z_right = c2_z + sign(slope) * sqrt(slope*slope/(1+slope*slope))*this->smoothing_radius;

  c1_r = this->rad_left + slope * ( z_left + this->length ) +
      sqrt( this->smoothing_radius * this->smoothing_radius  - SQR( z_left - c1_z ) );
  c2_r = this->rad_left + slope * ( z_right + this->length ) +
      sqrt( this->smoothing_radius * this->smoothing_radius  - SQR( z_right - c2_z ) );
  c1_r = this->rad_left+this->smoothing_radius;
  c2_r = this->rad_right+this->smoothing_radius;

  double c1_or = this->outer_rad_left-this->smoothing_radius;
  double c2_or = this->outer_rad_right-this->smoothing_radius;
 
  /* Check if we are in the region of the left wall */
  if (( (r >= c1_r) && (r <= c1_or) && (z <= c1_z) )) {
    dist_vector_z=-z - this->length;
    dist_vector_r=0;
    *dist = -z - this->length;
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return 0;
  }
  /* Check if we are in the region of the right wall */
  if (( (r >= c2_r) && (r<c2_or) && (z >= c2_z) ) ) {
    dist_vector_z=-z + this->length;
    dist_vector_r=0;
    *dist = +z - this->length;
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
    *dist=temp-this->smoothing_radius;
    dist_vector_r-=dist_vector_r/temp*this->smoothing_radius;
    dist_vector_z-=dist_vector_z/temp*this->smoothing_radius;
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return 0;
  }


  if ( p2_z>=c1_z && p2_z<=c2_z && dist_vector_r_o <= 0 ) {
    temp=sqrt( dist_vector_r_o*dist_vector_r_o + dist_vector_z_o*dist_vector_z_o );
    *dist=temp-this->smoothing_radius;
    dist_vector_r_o-=dist_vector_r_o/temp*this->smoothing_radius;
    dist_vector_z_o-=dist_vector_z_o/temp*this->smoothing_radius;
    for (i=0; i<3; i++) vec[i]=-dist_vector_r_o*e_r[i] - dist_vector_z_o*e_z[i];
    return 0;
  }


  /* Check if we are in the range of the left smoothing circle */
  if (p1_z <= c1_z && r <= c1_r ) {
    /* distance from the smoothing center */
    norm = sqrt( (z - c1_z)*(z - c1_z) + (r - c1_r)*(r - c1_r) );
    *dist = norm - this->smoothing_radius;
    dist_vector_r=(this->smoothing_radius/norm -1)*(r - c1_r);
    dist_vector_z=(this->smoothing_radius/norm - 1)*(z - c1_z);
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return 0;
  }
  /* upper left smoothing circle */
  if (p2_z <= c1_z && r >= c1_or ) {
    /* distance from the smoothing center */
    norm = sqrt( (z - c1_z)*(z - c1_z) + (r - c1_or)*(r - c1_or) );
    *dist = norm - this->smoothing_radius;
    dist_vector_r=(this->smoothing_radius/norm -1)*(r - c1_or);
    dist_vector_z=(this->smoothing_radius/norm - 1)*(z - c1_z);
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return  0;
  }
  /* Check if we are in the range of the right smoothing circle */
  if (p1_z >= c2_z && r <= c2_r ) {
    norm = sqrt( (z - c2_z)*(z - c2_z) + (r - c2_r)*(r - c2_r) );
    *dist = norm - this->smoothing_radius;
    dist_vector_r=(this->smoothing_radius/norm -1)*(r - c2_or);
    dist_vector_z=(this->smoothing_radius/norm - 1)*(z - c2_z);
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return  0;
  }
  /* Check if we are in the range of the upper right smoothing circle */
  if (p2_z >= c2_z && r >= c2_or ) {
    norm = sqrt( (z - c2_z)*(z - c2_z) + (r - c2_or)*(r - c2_or) );
    *dist = norm - this->smoothing_radius;
    dist_vector_r=(this->smoothing_radius/norm -1)*(r - c2_or);
    dist_vector_z=(this->smoothing_radius/norm - 1)*(z - c2_z);
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return  0;
  }
  *dist=-1e99;
  vec[0] = vec[1] = vec[2] = 1e99;
  return 0;
//  exit(printf("should never be reached, z %f, r%f\n",z, r));
}

int Maze::calculate_dist(double *ppos, double *dist, double *vec)
{
  int i,min_axis,cursph[3],dim;
  double diasph,fac,c_dist,sph_dist,cyl_dist,temp_dis;
  double sph_vec[3],cyl_vec[3];

  dim=(int) this->dim;
  diasph = box_l[0]/this->nsphere;

  /* First determine the distance to the sphere */
  c_dist=0.0;
  for(i=0;i<3;i++) {
    cursph[i] = (int) (ppos[i]/diasph);
    sph_vec[i] = (cursph[i]+0.5) * diasph  - (ppos[i]);
    c_dist += SQR(sph_vec[i]);
  }
  c_dist = sqrt(c_dist);
  sph_dist = this->sphrad - c_dist;
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
  cyl_dist = this->cylrad - c_dist;
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

int Cylinder::calculate_dist(double *ppos, double *dist, double *vec)
{
  int i;
  double d_per,d_par,d_real,d_per_vec[3],d_par_vec[3],d_real_vec[3];

  d_real = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = ppos[i] - this->pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);
    
  d_par=0.;
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * this->axis[i]);
  }
    
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * this->axis[i] ;
    d_per_vec[i] = ppos[i] - (this->pos[i] + d_par_vec[i]) ;
  }
		
  d_per=sqrt(SQR(d_real)-SQR(d_par));
  d_par = fabs(d_par) ;

  if ( this->direction == -1 ) {
    /*apply force towards inside cylinder */
    d_per = this->rad - d_per ;
    d_par = this->length - d_par;
    if (d_per < d_par )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= -d_per_vec[i] * d_per /  (this->rad - d_per) ;
      }
    } else {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= -d_par_vec[i] * d_par /  (this->length - d_par) ;
      }
    }
  } else {
    /*apply force towards outside cylinder */
    d_per = d_per - this->rad ;
    d_par = d_par - this->length ;
    if (d_par < 0 )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= d_per_vec[i] * d_per /  (d_per + this->rad) ;
      }
    } else if ( d_per < 0) {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= d_par_vec[i] * d_par /  (d_par + this->length) ;
      }
    } else {
      *dist = sqrt( SQR(d_par) + SQR(d_per)) ;
      for (i=0; i<3;i++) {
	vec[i]=
	  d_per_vec[i] * d_per /  (d_per + this->rad) +
	  d_par_vec[i] * d_par /  (d_par + this->length) ;
      }	
    }
  }
  return 0;
}

int SpheroCylinder::calculate_dist(double *ppos, double *dist, double *vec)
{
  int i;
  double d = 0.0;
  double ppos_local[3];

  for(i = 0; i < 3; i++) {
    ppos_local[i] = ppos[i] - this->pos[i];
    d += ppos_local[i] * this->axis[i];
  }

  if(abs(d) >= this->length) {
    *dist = 0.0;
    
    for(i = 0; i < 3; i++) {
      vec[i] = ppos_local[i] - this->length * this->axis[i] * sign(d);
      *dist += vec[i]*vec[i];
    }

    *dist = sqrt(*dist);
    
    if(*dist != 0.0)
      for(i = 0; i < 3; i++)
        vec[i] /= *dist;
    
    *dist -= this->rad;

    for(i = 0; i < 3; i++)
      vec[i] *= *dist;

    *dist *= this->direction;
  }
  else
    Cylinder::calculate_dist(ppos, dist, vec);
  return 0;
}

int Plane::calculate_dist(double* ppos, double *dist, double *vec)
{
  int i;
  double c_dist_sqr,c_dist;
  
  c_dist_sqr=0.0;
  for(i=0;i<3;i++) {
    if(this->pos[i] >= 0) {
      vec[i] = this->pos[i] - ppos[i];
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

