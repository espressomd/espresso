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
/** \file constraint.cpp
    Implementation of \ref constraint.hpp "constraint.hpp", here it's just the parsing stuff.
*/

#include <algorithm>
#include "constraint.hpp"
#include "energy_inline.hpp"
#include "forces_inline.hpp"
#include "tunable_slip.hpp"

// for the charged rod "constraint"
#define C_GAMMA   0.57721566490153286060651209008

int reflection_happened;

#ifdef CONSTRAINTS

int Constraint::n_constraints       = 0;
Constraint *Constraint::constraints = NULL;

Constraint *generate_constraint()
{
  Constraint::n_constraints++;
  Constraint::constraints = (Constraint*)realloc(Constraint::constraints,Constraint::n_constraints*sizeof(Constraint));
  Constraint::constraints[Constraint::n_constraints-1].type = CONSTRAINT_NONE;
  Constraint::constraints[Constraint::n_constraints-1].part_rep.p.identity = -Constraint::n_constraints;
  Constraint::constraints[Constraint::n_constraints-1].constraint_number = Constraint::n_constraints-1;
  
  return &Constraint::constraints[Constraint::n_constraints-1];
}

void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < Constraint::n_constraints; n++)
    for (i = 0; i < 3; i++)
      Constraint::constraints[n].part_rep.f.f[i] = 0;
}


static double sign(double x) {
  if (x > 0)
    return 1.;
  else
    return -1;
}

int Shape::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
{
  return 1; //if calculate_dist not implemented in derived class, calculate_dist fails
}

int Constraint_box::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
{
  *dist = -1;
  return 0;
}

int Constraint_wall::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
{
  int i;

  *dist = -this->d;
  for(i=0;i<3;i++) *dist += ppos[i]*this->n[i];
  
  for(i=0;i<3;i++) vec[i] = this->n[i] * *dist;
  return 0;
}


int Constraint_sphere::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
{
  int i;
  double fac,  c_dist;

  c_dist=0.0;
  for(i=0;i<3;i++) {
    vec[i] = this->pos[i] - ppos[i];
    c_dist += SQR(vec[i]);
  }
  c_dist = sqrt(c_dist);
  
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


int Constraint_maze::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
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

int Constraint_cylinder::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
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

int Constraint_spherocylinder::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
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
	  dynamic_cast<Constraint_cylinder*> (this)->calculate_dist(p1, ppos, c_p, dist, vec);
  return 0;
}

int Constraint_rhomboid::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec)
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

int Constraint_pore::calculate_dist(Particle *p1, double* ppos, Particle *c_p, double *dist, double *vec)
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

int Constraint_plane::calculate_dist(Particle *p1, double* ppos, Particle *c_p, double *dist, double *vec)
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

int Constraint_stomatocyte::calculate_dist( Particle *p1, double* ppos,
                                 Particle *c_p,
                                 double *dist, double *vec )
{
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

  a = this->outer_radius;
  b = this->inner_radius;
  c = this->layer_width;
  a = a*c;
  b = b*c;

  // Set the position and orientation of the stomatocyte 

  double stomatocyte_3D_position [3] = { this->position_x,
                                         this->position_y,
                                         this->position_z };

  double stomatocyte_3D_orientation [3] = { this->orientation_x,
                                            this->orientation_y,
                                            this->orientation_z };

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

  if ( this->direction == -1 )
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

int Constraint_hollow_cone::calculate_dist( Particle *p1, double ppos [3],
                                 Particle *c_p,
                                 double *dist, double *vec )
{
  // Parameters

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

  r0 = this->inner_radius;
  r1 = this->outer_radius;
  w = this->width;
  alpha = this->opening_angle;

  // Set the position and orientation of the hollow cone

  double hollow_cone_3D_position [3] = { this->position_x,
                                         this->position_y,
                                         this->position_z };

  double hollow_cone_3D_orientation [3] = { this->orientation_x,
                                            this->orientation_y,
                                            this->orientation_z };

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

  if ( this->direction == -1 )
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

void Constraint_rod::add_constraint_force_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* force, double* torque1, double* torque2, IA_parameters* ia_params)
{
#ifdef ELECTROSTATICS
  int i;
  double fac, vec[2], c_dist_2;
  Particle* c_p = &current_constraint->part_rep;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = folded_pos[i] - this->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && this->lambda != 0.0) {
    fac = 2*coulomb.prefactor*this->lambda*p1->p.q/c_dist_2;
    p1->f.f[0]  += fac*vec[0];
    p1->f.f[1]  += fac*vec[1];
    c_p->f.f[0] -= fac*vec[0];
    c_p->f.f[1] -= fac*vec[1];
  }
#endif
}

void Constraint_rod::add_constraint_energy_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* nonbonded_en, double* coulomb_en, double* magnetic_en, IA_parameters* ia_params)
{
#ifdef ELECTROSTATICS
  int i;
  double vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = folded_pos[i] - this->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && this->lambda != 0.0)
    *coulomb_en = coulomb.prefactor*p1->p.q*this->lambda*(-log(c_dist_2*SQR(box_l_i[2])) + 2*(M_LN2 - C_GAMMA));
#endif
}

void Constraint_plate::add_constraint_force_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* force, double* torque1, double* torque2, IA_parameters* ia_params)
{
#ifdef ELECTROSTATICS

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && this->sigma != 0.0) {
    *force = 2*M_PI*coulomb.prefactor*this->sigma*p1->p.q;
    if (folded_pos[2] < this->pos)
      *force = -*force;
  }
#endif
}

void Constraint_plate::add_constraint_energy_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* nonbonded_en, double* coulomb_en, double* magnetic_en, IA_parameters* ia_params)
{
#ifdef ELECTROSTATICS
  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && this->sigma != 0.0)
    *coulomb_en = -2*M_PI*coulomb.prefactor*this->sigma*p1->p.q*fabs(folded_pos[2] - this->pos);
#endif

}

//ER
void Constraint_ext_magn_field::add_constraint_force_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* force, double* torque1, double* torque2, IA_parameters* ia_params)
{
#ifdef ROTATION
#ifdef DIPOLES
  p1->f.torque[0] += p1->r.dip[1]*this->ext_magn_field[2]-p1->r.dip[2]*this->ext_magn_field[1];
  p1->f.torque[1] += p1->r.dip[2]*this->ext_magn_field[0]-p1->r.dip[0]*this->ext_magn_field[2];
  p1->f.torque[2] += p1->r.dip[0]*this->ext_magn_field[1]-p1->r.dip[1]*this->ext_magn_field[0];
#endif
#endif
}

void Constraint_ext_magn_field::add_constraint_energy_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* nonbonded_en, double* coulomb_en, double* magnetic_en, IA_parameters* ia_params)
{
#ifdef DIPOLES
     *magnetic_en = -1.0 * scalar(this->ext_magn_field,p1->r.dip);
#endif
}
//end ER

void reflect_particle(Particle *p1, double *distance_vec, int reflecting) {
  double vec[3];
  double norm; 

      double folded_pos[3];
      int img[3];

      memcpy(folded_pos, p1->r.p, 3*sizeof(double));
      memcpy(img, p1->l.i, 3*sizeof(int));
      fold_position(folded_pos, img);

      memcpy(vec, distance_vec, 3*sizeof(double));
/* For Debugging your can show the folded coordinates of the particle before
 * and after the reflecting by uncommenting these lines  */
 //     printf("position before reflection %f %f %f\n",folded_pos[0], folded_pos[1], folded_pos[2]); 


       reflection_happened = 1;
       norm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
       p1->r.p[0] = p1->r.p[0]-2*vec[0];
       p1->r.p[1] = p1->r.p[1]-2*vec[1];
       p1->r.p[2] = p1->r.p[2]-2*vec[2];

   /*  This can show the folded position after reflection      
       memcpy(folded_pos, p1->r.p, 3*sizeof(double));
       memcpy(img, p1->l.i, 3*sizeof(int));
       fold_position(folded_pos, img);
       printf("position after reflection %f %f %f\n",folded_pos[0], folded_pos[1], folded_pos[2]); */

       /* vec seams to be the vector that points from the wall to the particle*/
       /* now normalize it */ 
       if ( reflecting==1 ) {
         vec[0] /= norm;
         vec[1] /= norm;
         vec[2] /= norm;
         /* calculating scalar product - reusing var norm */
         norm = vec[0] *  p1->m.v[0] + vec[1] * p1->m.v[1] + vec[2] * p1->m.v[2];
         /* now add twice the normal component to the velcity */
          p1->m.v[0] = p1->m.v[0]-2*vec[0]*norm; /* norm is still the scalar product! */
          p1->m.v[1] = p1->m.v[1]-2*vec[1]*norm;
          p1->m.v[2] = p1->m.v[2]-2*vec[2]*norm;
       } else if (reflecting == 2) {
         /* if bounce back, invert velocity */
          p1->m.v[0] =-p1->m.v[0]; /* norm is still the scalar product! */
          p1->m.v[1] =-p1->m.v[1];
          p1->m.v[2] =-p1->m.v[2];

       }

}

void Shape::add_constraint_force_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* force, double* torque1, double* torque2, IA_parameters* ia_params)
{
	double dist, vec[3];
	if(checkIfInteraction(ia_params)) {
		if (!(current_constraint->_shape->calculate_dist(p1, folded_pos, &current_constraint->part_rep, &dist, vec))) {
			if ( dist > 0 ) {
				calc_non_bonded_pair_force(p1, &current_constraint->part_rep,
						ia_params,vec,dist,dist*dist, force,
						torque1, torque2);
#ifdef TUNABLE_SLIP
				if (current_constraint->_shape->tunable_slip == 1 ) {
					add_tunable_slip_pair_force(p1, &current_constraint->part_rep,ia_params,vec,dist,force);
				}
#endif
			}
			else if ( dist <= 0 && current_constraint->_shape->penetrable == 1 ) {
				if ( (current_constraint->_shape->only_positive != 1) && ( dist < 0 ) ) {
					calc_non_bonded_pair_force(p1, &current_constraint->part_rep,
							ia_params,vec,-1.0*dist,dist*dist, force,
							torque1, torque2);
				}
			}
			else {
				if(current_constraint->_shape->reflecting){
					reflect_particle(p1, &(vec[0]), current_constraint->_shape->reflecting);
				} else {
					ostringstream msg;
					msg <<"constraint "<< current_constraint->constraint_number<<" violated by particle "<<p1->p.identity;
					runtimeError(msg);
				}
			}
		}
		else {
			ostringstream msg;
			msg <<"ERROR: encountered unknown constraint during force computation for constraint"<< current_constraint->constraint_number;
			runtimeError(msg);

		}
	}
}

void add_constraints_forces(Particle *p1)
{
	if (Constraint::n_constraints==0)
		return;
	int n, j;
	double force[3], torque1[3], torque2[3];

	IA_parameters *ia_params;
	double folded_pos[3];
	int img[3];

	/* fold the coordinate[2] of the particle */
	memcpy(folded_pos, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(folded_pos, img);

	for(n=0;n<Constraint::n_constraints;n++) {
		ia_params=get_ia_param(p1->p.type, (&Constraint::constraints[n].part_rep)->p.type);
		for (j = 0; j < 3; j++) {
			force[j] = 0;
#ifdef ROTATION
			torque1[j] = torque2[j] = 0;
#endif
		}
		Constraint::constraints[n]._shape->add_constraint_force_default(&Constraint::constraints[n],p1,folded_pos,force,torque1,torque2, ia_params);

		for (j = 0; j < 3; j++) {
			p1->f.f[j] += force[j];
			Constraint::constraints[n].part_rep.f.f[j] -= force[j];
#ifdef ROTATION
			p1->f.torque[j] += torque1[j];
			Constraint::constraints[n].part_rep.f.torque[j] += torque2[j];
#endif
		}
	}
}
void Shape::add_constraint_energy_default (Constraint* current_constraint, Particle *p1, double* folded_pos, double* nonbonded_en, double* coulomb_en, double* magnetic_en, IA_parameters* ia_params)
{
	double dist, vec[3];
	magnetic_en = 0;
	coulomb_en = 0;
	if(checkIfInteraction(ia_params)) {
		if (!(current_constraint->_shape->calculate_dist(p1, folded_pos, &current_constraint->part_rep, &dist, vec))) {
			if ( dist > 0 ) {
				*nonbonded_en = calc_non_bonded_pair_energy(p1, &current_constraint->part_rep,
						ia_params,vec,dist,dist*dist);
			}
			else if ( dist <= 0 && current_constraint->_shape->penetrable == 1 ) {
				if ( (current_constraint->_shape->only_positive != 1) && ( dist < 0 ) ) {
					*nonbonded_en = calc_non_bonded_pair_energy(p1, &current_constraint->part_rep,
							ia_params,vec,-1.0*dist,dist*dist);
				}
			}
			else {
					ostringstream msg;
					msg <<"constraint "<< current_constraint->constraint_number<<" violated by particle making the energy not calculable"<<p1->p.identity;
					runtimeError(msg);
			}
		}
		else {
			fprintf (stderr,"Warning:encountered a constraint for which the energy cannot be computed\n");
		}
	}
}
double add_constraints_energy(Particle *p1)
{
	int n, type;
	double nonbonded_en, coulomb_en, magnetic_en;
	IA_parameters *ia_params;
	double folded_pos[3];
	int img[3];

	/* fold the coordinate[2] of the particle */
	memcpy(folded_pos, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(folded_pos, img);
	for(n=0;n<Constraint::n_constraints;n++) {
		ia_params = get_ia_param(p1->p.type, (&Constraint::constraints[n].part_rep)->p.type);
		nonbonded_en = 0.;
		coulomb_en   = 0.;
		magnetic_en = 0.;

		Constraint::constraints[n]._shape->add_constraint_energy_default (&Constraint::constraints[n], p1, folded_pos, &nonbonded_en, &coulomb_en, &magnetic_en, ia_params);
		if (energy.n_coulomb > 0)
			energy.coulomb[0] += coulomb_en;

		if (energy.n_dipolar > 0)
			energy.dipolar[0] += magnetic_en;

		type = (&Constraint::constraints[n].part_rep)->p.type;
		if (type >= 0)
			*obsstat_nonbonded(&energy, p1->p.type, type) += nonbonded_en;
	}
	return 0.;
}

int Constraint_slitpore::calculate_dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec) {
  // the left circles
  double box_l_x = box_l[0];
  double c11[2] = { box_l_x/2-this->pore_width/2-this->upper_smoothing_radius, this->pore_mouth - this->upper_smoothing_radius };
  double c12[2] = { box_l_x/2-this->pore_width/2+this->lower_smoothing_radius, this->pore_mouth - this->pore_length  + this->lower_smoothing_radius };
  // the right circles
  double c21[2] = { box_l_x/2+this->pore_width/2+this->upper_smoothing_radius, this->pore_mouth - this->upper_smoothing_radius };
  double c22[2] = { box_l_x/2+this->pore_width/2-this->lower_smoothing_radius, this->pore_mouth - this->pore_length  + this->lower_smoothing_radius };

//  printf("c11 %f %f\n", c11[0], c11[1]);
//  printf("c12 %f %f\n", c12[0], c12[1]);
//  printf("c21 %f %f\n", c21[0], c21[1]);
//  printf("c22 %f %f\n", c22[0], c22[1]);


  if (ppos[2] > this->pore_mouth + this->channel_width/2) {
//    printf("upper wall\n");
    // Feel the upper wall
    *dist = this->pore_mouth + this->channel_width - ppos[2];
    vec[0] = vec[1] = 0;
    vec[2] = -*dist;
    return 0;
  }

  if (ppos[0]<c11[0] || ppos[0] > c21[0]) {
    // Feel the lower wall of the channel
//    printf("lower wall\n");
    *dist = ppos[2] - this->pore_mouth;
    vec[0] = vec[1] = 0;
    vec[2] = *dist;
    return 0;
  }

  if (ppos[2] > c11[1]) {
    // Feel the upper smoothing
    if (ppos[0] < box_l_x/2) {
//    printf("upper smoothing left\n");
      *dist = sqrt( SQR(c11[0] - ppos[0]) + SQR(c11[1] - ppos[2])) - this->upper_smoothing_radius;
      vec[0] = -( c11[0] - ppos[0] ) * (*dist)/(*dist+this->upper_smoothing_radius);
      vec[1] = 0;
      vec[2] = -( c11[1] - ppos[2] ) * (*dist)/(*dist+this->upper_smoothing_radius);
      return 0;
    } else {
//    printf("upper smoothing right\n");
      *dist = sqrt( SQR(c21[0] - ppos[0]) + SQR(c21[1] - ppos[2])) - this->upper_smoothing_radius;
      vec[0] = -( c21[0] - ppos[0] ) * (*dist)/(*dist+this->upper_smoothing_radius);
      vec[1] = 0;
      vec[2] = -( c21[1] - ppos[2] ) * (*dist)/(*dist+this->upper_smoothing_radius);
      return 0;
    }
  }
  
  if (ppos[2] > c12[1]) {
    // Feel the pore wall
    if (ppos[0] < box_l_x/2) {
//    printf("pore left\n");
      *dist = ppos[0] - (box_l_x/2-this->pore_width/2);
      vec[0]=*dist;
      vec[1]=vec[2]=0;
      return 0;
    } else {
//    printf("pore right\n");
      *dist =  (box_l_x/2+this->pore_width/2) - ppos[0];
      vec[0]=-*dist;
      vec[1]=vec[2]=0;
      return 0;
    }
  }

  if (ppos[0]>c12[0] && ppos[0] < c22[0]) {
//    printf("pore end\n");
    // Feel the pore end wall
    *dist = ppos[2] - (this->pore_mouth-this->pore_length);
    vec[0]=vec[1]=0;
    vec[2]=*dist;
    return 0;
  }
  // Else
  // Feel the lower smoothing
    if (ppos[0] < box_l_x/2) {
//    printf("lower smoothing left\n");
      *dist = -sqrt( SQR(c12[0] - ppos[0]) + SQR(c12[1] - ppos[2])) + this->lower_smoothing_radius;
      vec[0] = ( c12[0] - ppos[0] ) * (*dist)/(-*dist+this->lower_smoothing_radius);
      vec[1] = 0;
      vec[2] = ( c12[1] - ppos[2] ) * (*dist)/(-*dist+this->lower_smoothing_radius);
      return 0;
    } else {
//    printf("lower smoothing right\n");
      *dist = -sqrt( SQR(c22[0] - ppos[0]) + SQR(c22[1] - ppos[2])) + this->lower_smoothing_radius;
      vec[0] = ( c22[0] - ppos[0] ) * (*dist)/(-*dist+this->lower_smoothing_radius);
      vec[1] = 0;
      vec[2] = ( c22[1] - ppos[2] ) * (*dist)/(-*dist+this->lower_smoothing_radius);
      return 0;
    }

    return 0;
}
#endif

