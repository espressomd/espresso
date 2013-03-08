/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file constraint.c
    Implementation of \ref constraint.h "constraint.h", here it's just the parsing stuff.
*/

#include "constraint.h"
#include "energy.h"
#include "forces.h"
#include "tunable_slip.h"

// for the charged rod "constraint"
#define C_GAMMA   0.57721566490153286060651209008

int reflection_happened;

#ifdef CONSTRAINTS

int n_constraints       = 0;
Constraint *constraints = NULL;

Constraint *generate_constraint()
{
  n_constraints++;
  constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  constraints[n_constraints-1].type = CONSTRAINT_NONE;
  constraints[n_constraints-1].part_rep.p.identity = -n_constraints;
  
  return &constraints[n_constraints-1];
}

void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < n_constraints; n++)
    for (i = 0; i < 3; i++)
      constraints[n].part_rep.f.f[i] = 0;
}


static double sign(double x) {
  if (x > 0)
    return 1.;
  else
    return -1;
}

static double max(double x1, double x2) {
  return x1>x2?x1:x2;
}

void calculate_wall_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_wall *c, double *dist, double *vec)
{
  int i;

  *dist = -c->d;
  for(i=0;i<3;i++) *dist += ppos[i]*c->n[i];
  
  for(i=0;i<3;i++) vec[i] = c->n[i] * *dist;
  
}


void calculate_sphere_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_sphere *c, double *dist, double *vec)
{
  int i;
  double fac,  c_dist;

  c_dist=0.0;
  for(i=0;i<3;i++) {
    vec[i] = c->pos[i] - ppos[i];
    c_dist += SQR(vec[i]);
  }
  c_dist = sqrt(c_dist);
  
  if ( c->direction == -1 ) {
  /* apply force towards inside the sphere */
    *dist = c->rad - c_dist;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= fac;
  } else {
   /* apply force towards outside the sphere */
    *dist = c_dist - c->rad;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= -fac;
  }
}


void calculate_maze_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_maze *c, double *dist, double *vec)
{
  int i,min_axis,cursph[3],dim;
  double diasph,fac,c_dist,sph_dist,cyl_dist,temp_dis;
  double sph_vec[3],cyl_vec[3];

  dim=(int) c->dim;
  diasph = box_l[0]/c->nsphere;

  /* First determine the distance to the sphere */
  c_dist=0.0;
  for(i=0;i<3;i++) {
    cursph[i] = (int) (ppos[i]/diasph);
    sph_vec[i] = (cursph[i]+0.5) * diasph  - (ppos[i]);
    c_dist += SQR(sph_vec[i]);
  }
  c_dist = sqrt(c_dist);
  sph_dist = c->sphrad - c_dist;
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
  cyl_dist = c->cylrad - c_dist;
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
}

void calculate_cylinder_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_cylinder *c, double *dist, double *vec)
{
  int i;
  double d_per,d_par,d_real,d_per_vec[3],d_par_vec[3],d_real_vec[3];

  d_real = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = ppos[i] - c->pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);
    
  d_par=0.;
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * c->axis[i]);
  }
    
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * c->axis[i] ;
    d_per_vec[i] = ppos[i] - (c->pos[i] + d_par_vec[i]) ;
  }
		
  d_per=sqrt(SQR(d_real)-SQR(d_par));
  d_par = fabs(d_par) ;

  if ( c->direction == -1 ) {
    /*apply force towards inside cylinder */
    d_per = c->rad - d_per ;
    d_par = c->length - d_par;
    if (d_per < d_par )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= -d_per_vec[i] * d_per /  (c->rad - d_per) ;
      }
    } else {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= -d_par_vec[i] * d_par /  (c->length - d_par) ;
      }
    }
  } else {
    /*apply force towards outside cylinder */
    d_per = d_per - c->rad ;
    d_par = d_par - c->length ;
    if (d_par < 0 )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= d_per_vec[i] * d_per /  (d_per + c->rad) ;
      }
    } else if ( d_per < 0) {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= d_par_vec[i] * d_par /  (d_par + c->length) ;
      }
    } else {
      *dist = sqrt( SQR(d_par) + SQR(d_per)) ;
      for (i=0; i<3;i++) {
	vec[i]=
	  d_per_vec[i] * d_per /  (d_per + c->rad) +
	  d_par_vec[i] * d_par /  (d_par + c->length) ;
      }	
    }
  }
}

void calculate_rhomboid_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_rhomboid *c, double *dist, double *vec)
{	
	double axb[3], bxc[3], axc[3];
	double A, B, C;
	double a_dot_bxc, b_dot_axc, c_dot_axb;
	double tmp;
	double d;
	
	//calculate a couple of vectors and scalars that are going to be used frequently
	
	axb[0] = c->a[1]*c->b[2] - c->a[2]*c->b[1];
	axb[1] = c->a[2]*c->b[0] - c->a[0]*c->b[2];
	axb[2] = c->a[0]*c->b[1] - c->a[1]*c->b[0];
	
	bxc[0] = c->b[1]*c->c[2] - c->b[2]*c->c[1];
	bxc[1] = c->b[2]*c->c[0] - c->b[0]*c->c[2];
	bxc[2] = c->b[0]*c->c[1] - c->b[1]*c->c[0];
	
	axc[0] = c->a[1]*c->c[2] - c->a[2]*c->c[1];
	axc[1] = c->a[2]*c->c[0] - c->a[0]*c->c[2];
	axc[2] = c->a[0]*c->c[1] - c->a[1]*c->c[0];
	
	a_dot_bxc = c->a[0]*bxc[0] + c->a[1]*bxc[1] + c->a[2]*bxc[2];
	b_dot_axc = c->b[0]*axc[0] + c->b[1]*axc[1] + c->b[2]*axc[2];
	c_dot_axb = c->c[0]*axb[0] + c->c[1]*axb[1] + c->c[2]*axb[2];
	
	//represent the distance from pos to ppos as a linear combination of the edge vectors.
	
	A = (ppos[0]-c->pos[0])*bxc[0] + (ppos[1]-c->pos[1])*bxc[1] + (ppos[2]-c->pos[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0])*axc[0] + (ppos[1]-c->pos[1])*axc[1] + (ppos[2]-c->pos[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0])*axb[0] + (ppos[1]-c->pos[1])*axb[1] + (ppos[2]-c->pos[2])*axb[2];
	C /= c_dot_axb;
	
	//the coefficients tell whether ppos lies within the cone defined by pos and the adjacent edges
	
	if(A <= 0 && B <= 0 && C <= 0)
	{
		vec[0] = ppos[0]-c->pos[0];
		vec[1] = ppos[1]-c->pos[1];
		vec[2] = ppos[2]-c->pos[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
	  return;
	}
	
	//check for cone at pos+a

	A = (ppos[0]-c->pos[0]-c->a[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->a[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->a[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2])*axb[2];
	C /= c_dot_axb;

	if(A >= 0 && B <= 0 && C <= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->a[0];
		vec[1] = ppos[1]-c->pos[1]-c->a[1];
		vec[2] = ppos[2]-c->pos[2]-c->a[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
	  return;
	}
	
	//check for cone at pos+b

	A = (ppos[0]-c->pos[0]-c->b[0])*bxc[0] + (ppos[1]-c->pos[1]-c->b[1])*bxc[1] + (ppos[2]-c->pos[2]-c->b[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->b[0])*axc[0] + (ppos[1]-c->pos[1]-c->b[1])*axc[1] + (ppos[2]-c->pos[2]-c->b[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->b[0])*axb[0] + (ppos[1]-c->pos[1]-c->b[1])*axb[1] + (ppos[2]-c->pos[2]-c->b[2])*axb[2];
	C /= c_dot_axb;

	if(A <= 0 && B >= 0 && C <= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->b[0];
		vec[1] = ppos[1]-c->pos[1]-c->b[1];
		vec[2] = ppos[2]-c->pos[2]-c->b[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
		return;
	}
	
	//check for cone at pos+c

	A = (ppos[0]-c->pos[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->c[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(A <= 0 && B <= 0 && C >= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->c[0];
		vec[1] = ppos[1]-c->pos[1]-c->c[1];
		vec[2] = ppos[2]-c->pos[2]-c->c[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	
	  return;
	}
	
	//check for cone at pos+a+b

	A = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*axb[2];
	C /= c_dot_axb;

	if(A >= 0 && B >= 0 && C <= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->b[0];
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->b[1];
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->b[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for cone at pos+a+c

	A = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(A >= 0 && B <= 0 && C >= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->c[0];
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->c[1];
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->c[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for cone at pos+a+c

	A = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(A <= 0 && B >= 0 && C >= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->b[0]-c->c[0];
		vec[1] = ppos[1]-c->pos[1]-c->b[1]-c->c[1];
		vec[2] = ppos[2]-c->pos[2]-c->b[2]-c->c[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for cone at pos+a+b+c

	A = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;	
	B = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*axc[2];
	B /= b_dot_axc;	
	C = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(A >= 0 && B >= 0 && C >= 0)
	{
		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0];
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1];
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2];
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos, a
	
	B = (ppos[0]-c->pos[0])*axc[0] + (ppos[1]-c->pos[1])*axc[1] + (ppos[2]-c->pos[2])*axc[2];
	B /= b_dot_axc;
	C = (ppos[0]-c->pos[0])*axb[0] + (ppos[1]-c->pos[1])*axb[1] + (ppos[2]-c->pos[2])*axb[2];
	C /= c_dot_axb;
	
	if(B <= 0 && C <= 0)
	{
		tmp = (ppos[0]-c->pos[0])*c->a[0] + (ppos[1]-c->pos[1])*c->a[1] + (ppos[2]-c->pos[2])*c->a[2];
		tmp /= c->a[0]*c->a[0] + c->a[1]*c->a[1] + c->a[2]*c->a[2];
		
		vec[0] = ppos[0]-c->pos[0] - c->a[0]*tmp;
		vec[1] = ppos[1]-c->pos[1] - c->a[1]*tmp;
		vec[2] = ppos[2]-c->pos[2] - c->a[2]*tmp;
		
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		
    return;
	}
	
	//check for prism at edge pos, b

	A = (ppos[0]-c->pos[0])*bxc[0] + (ppos[1]-c->pos[1])*bxc[1] + (ppos[2]-c->pos[2])*bxc[2];
	A /= a_dot_bxc;
	C = (ppos[0]-c->pos[0])*axb[0] + (ppos[1]-c->pos[1])*axb[1] + (ppos[2]-c->pos[2])*axb[2];
	C /= c_dot_axb;

	if(A <= 0 && C <= 0)
	{
		tmp = (ppos[0]-c->pos[0])*c->b[0] + (ppos[1]-c->pos[1])*c->b[1] + (ppos[2]-c->pos[2])*c->b[2];
		tmp /= c->b[0]*c->b[0] + c->b[1]*c->b[1] + c->b[2]*c->b[2];
	
		vec[0] = ppos[0]-c->pos[0] - c->b[0]*tmp;
		vec[1] = ppos[1]-c->pos[1] - c->b[1]*tmp;
		vec[2] = ppos[2]-c->pos[2] - c->b[2]*tmp;
	
		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	
    return;
	}
	
	//check for prism at edge pos, c

	A = (ppos[0]-c->pos[0])*bxc[0] + (ppos[1]-c->pos[1])*bxc[1] + (ppos[2]-c->pos[2])*bxc[2];
	A /= a_dot_bxc;
	B = (ppos[0]-c->pos[0])*axc[0] + (ppos[1]-c->pos[1])*axc[1] + (ppos[2]-c->pos[2])*axc[2];
	B /= b_dot_axc;

	if(A <= 0 && B <= 0)
	{
		tmp = (ppos[0]-c->pos[0])*c->c[0] + (ppos[1]-c->pos[1])*c->c[1] + (ppos[2]-c->pos[2])*c->c[2];
		tmp /= c->c[0]*c->c[0] + c->c[1]*c->c[1] + c->c[2]*c->c[2];

		vec[0] = ppos[0]-c->pos[0] - c->c[0]*tmp;
		vec[1] = ppos[1]-c->pos[1] - c->c[1]*tmp;
		vec[2] = ppos[2]-c->pos[2] - c->c[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+a, b

	A = (ppos[0]-c->pos[0]-c->a[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2])*bxc[2];
	A /= a_dot_bxc;
	C = (ppos[0]-c->pos[0]-c->a[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2])*axb[2];
	C /= c_dot_axb;

	if(A >= 0 && C <= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->a[0])*c->b[0] + (ppos[1]-c->pos[1]-c->a[1])*c->b[1] + (ppos[2]-c->pos[2]-c->a[2])*c->b[2];
		tmp /= c->b[0]*c->b[0] + c->b[1]*c->b[1] + c->b[2]*c->b[2];

		vec[0] = ppos[0]-c->pos[0]-c->a[0] - c->b[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->a[1] - c->b[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->a[2] - c->b[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+a, c

	A = (ppos[0]-c->pos[0]-c->a[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2])*bxc[2];
	A /= a_dot_bxc;
	B = (ppos[0]-c->pos[0]-c->a[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2])*axc[2];
	B /= b_dot_axc;

	if(A >= 0 && B <= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->a[0])*c->c[0] + (ppos[1]-c->pos[1]-c->a[1])*c->c[1] + (ppos[2]-c->pos[2]-c->a[2])*c->c[2];
		tmp /= c->c[0]*c->c[0] + c->c[1]*c->c[1] + c->c[2]*c->c[2];

		vec[0] = ppos[0]-c->pos[0]-c->a[0] - c->c[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->a[1] - c->c[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->a[2] - c->c[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+b+c, c

	A = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;
	B = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*axc[2];
	B /= b_dot_axc;

	if(A <= 0 && B >= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*c->c[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*c->c[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*c->c[2];
		tmp /= c->c[0]*c->c[0] + c->c[1]*c->c[1] + c->c[2]*c->c[2];

		vec[0] = ppos[0]-c->pos[0]-c->b[0]-c->c[0] - c->c[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->b[1]-c->c[1] - c->c[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->b[2]-c->c[2] - c->c[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+b+c, b

	A = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;
	C = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(A <= 0 && C >= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*c->b[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*c->b[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*c->b[2];
		tmp /= c->b[0]*c->b[0] + c->b[1]*c->b[1] + c->b[2]*c->b[2];

		vec[0] = ppos[0]-c->pos[0]-c->b[0]-c->c[0] - c->b[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->b[1]-c->c[1] - c->b[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->b[2]-c->c[2] - c->b[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+b+c, a

	B = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*axc[2];
	B /= b_dot_axc;
	C = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*axb[2];
	C /= c_dot_axb;
																
	if(B >= 0 && C >= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->b[0]-c->c[0])*c->a[0] + (ppos[1]-c->pos[1]-c->b[1]-c->c[1])*c->a[1] + (ppos[2]-c->pos[2]-c->b[2]-c->c[2])*c->a[2];
		tmp /= c->a[0]*c->a[0] + c->a[1]*c->a[1] + c->a[2]*c->a[2];

		vec[0] = ppos[0]-c->pos[0]-c->b[0]-c->c[0] - c->a[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->b[1]-c->c[1] - c->a[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->b[2]-c->c[2] - c->a[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+a+b, a

	B = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*axc[2];
	B /= b_dot_axc;
	C = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*axb[2];
	C /= c_dot_axb;

	if(B >= 0 && C <= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*c->a[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*c->a[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*c->a[2];
		tmp /= c->a[0]*c->a[0] + c->a[1]*c->a[1] + c->a[2]*c->a[2];

		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->b[0] - c->a[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->b[1] - c->a[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->b[2] - c->a[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+a+b, c

	A = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*bxc[2];
	A /= a_dot_bxc;
	B = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*axc[2];
	B /= b_dot_axc;

	if(A >= 0 && B >= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->a[0]-c->b[0])*c->c[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1])*c->c[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2])*c->c[2];
		tmp /= c->c[0]*c->c[0] + c->c[1]*c->c[1] + c->c[2]*c->c[2];

		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->b[0] - c->c[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->b[1] - c->c[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->b[2] - c->c[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+a+c, a

	B = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*axc[2];
	B /= b_dot_axc;
	C = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(B <= 0 && C >= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*c->a[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*c->a[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*c->a[2];
		tmp /= c->a[0]*c->a[0] + c->a[1]*c->a[1] + c->a[2]*c->a[2];

		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->c[0] - c->a[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->c[1] - c->a[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->c[2] - c->a[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for prism at edge pos+a+c, b

	A = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*bxc[2];
	A /= a_dot_bxc;
	C = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*axb[2];
	C /= c_dot_axb;

	if(A >= 0 && C >= 0)
	{
		tmp = (ppos[0]-c->pos[0]-c->a[0]-c->c[0])*c->b[0] + (ppos[1]-c->pos[1]-c->a[1]-c->c[1])*c->b[1] + (ppos[2]-c->pos[2]-c->a[2]-c->c[2])*c->b[2];
		tmp /= c->b[0]*c->b[0] + c->b[1]*c->b[1] + c->b[2]*c->b[2];

		vec[0] = ppos[0]-c->pos[0]-c->a[0]-c->c[0] - c->b[0]*tmp;
		vec[1] = ppos[1]-c->pos[1]-c->a[1]-c->c[1] - c->b[1]*tmp;
		vec[2] = ppos[2]-c->pos[2]-c->a[2]-c->c[2] - c->b[2]*tmp;

		*dist = c->direction * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    return;
	}
	
	//check for face with normal -axb
	
	*dist = (ppos[0]-c->pos[0])*axb[0] + (ppos[1]-c->pos[1])*axb[1] + (ppos[2]-c->pos[2])*axb[2];
	*dist *= -1.;
	
	if(*dist >= 0)
	{
		tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
		*dist /= tmp;
	
		vec[0] = -*dist * axb[0]/tmp;
		vec[1] = -*dist * axb[1]/tmp;
		vec[2] = -*dist * axb[2]/tmp;
	
		*dist *= c->direction;

    return;
	}
	
	//calculate distance to face with normal axc

	*dist = (ppos[0]-c->pos[0])*axc[0] + (ppos[1]-c->pos[1])*axc[1] + (ppos[2]-c->pos[2])*axc[2];
	
	if(*dist >= 0)
	{
		tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
		*dist /= tmp;
	
		vec[0] = *dist * axc[0]/tmp;
		vec[1] = *dist * axc[1]/tmp;
		vec[2] = *dist * axc[2]/tmp;

		*dist *= c->direction;

    return;
	}
	
	//calculate distance to face with normal -bxc

	*dist = (ppos[0]-c->pos[0])*bxc[0] + (ppos[1]-c->pos[1])*bxc[1] + (ppos[2]-c->pos[2])*bxc[2];
	*dist *= -1.;
	
	if(*dist >= 0)
	{
		tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
		*dist /= tmp;
		
		vec[0] = -*dist * bxc[0]/tmp;
		vec[1] = -*dist * bxc[1]/tmp;
		vec[2] = -*dist * bxc[2]/tmp;
		
		*dist *= c->direction;

    return;
	}
	
	//calculate distance to face with normal axb
	
	*dist = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*axb[2];
	
	if(*dist >= 0)
	{
		tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
		*dist /= tmp;

		vec[0] = *dist * axb[0]/tmp;
		vec[1] = *dist * axb[1]/tmp;
		vec[2] = *dist * axb[2]/tmp;
	
		*dist *= c->direction;

    return;
	}
																					
	//calculate distance to face with normal -axc

	*dist = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*axc[2];
	*dist *= -1.;
	
	if(*dist >= 0)
	{
		tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
		*dist /= tmp;

		vec[0] = -*dist * axc[0]/tmp;
		vec[1] = -*dist * axc[1]/tmp;
		vec[2] = -*dist * axc[2]/tmp;
		
		*dist *= c->direction;

    return;
	}
																		
	//calculate distance to face with normal bxc

	*dist = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*bxc[2];
	
	if(*dist >= 0)
	{
		tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
		*dist /= tmp;

		vec[0] = *dist * bxc[0]/tmp;
		vec[1] = *dist * bxc[1]/tmp;
		vec[2] = *dist * bxc[2]/tmp;
		
		*dist *= c->direction;

    return;
	}
	
	//ppos lies within rhomboid. Find nearest wall for interaction.
	 
	//check for face with normal -axb
	
	*dist = (ppos[0]-c->pos[0])*axb[0] + (ppos[1]-c->pos[1])*axb[1] + (ppos[2]-c->pos[2])*axb[2];
	*dist *= -1.;
	tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
	*dist /= tmp;
	
	vec[0] = -*dist * axb[0]/tmp;
	vec[1] = -*dist * axb[1]/tmp;
	vec[2] = -*dist * axb[2]/tmp;
	
	*dist *= c->direction;

	//calculate distance to face with normal axc

	d = (ppos[0]-c->pos[0])*axc[0] + (ppos[1]-c->pos[1])*axc[1] + (ppos[2]-c->pos[2])*axc[2];
	tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
	d /= tmp;
	
	if(abs(d) < abs(*dist))
	{
		vec[0] = d * axc[0]/tmp;
		vec[1] = d * axc[1]/tmp;
		vec[2] = d * axc[2]/tmp;
	
		*dist = c->direction * d;
	}

	//calculate distance to face with normal -bxc

	d = (ppos[0]-c->pos[0])*bxc[0] + (ppos[1]-c->pos[1])*bxc[1] + (ppos[2]-c->pos[2])*bxc[2];
	d *= -1.;
	tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
	d /= tmp;

	if(abs(d) < abs(*dist))
	{							
		vec[0] = -d * bxc[0]/tmp;
		vec[1] = -d * bxc[1]/tmp;
		vec[2] = -d * bxc[2]/tmp;

		*dist = c->direction * d;
	}
	
	//calculate distance to face with normal axb

	d = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*axb[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*axb[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*axb[2];
	tmp = sqrt( axb[0]*axb[0] + axb[1]*axb[1] + axb[2]*axb[2] );
	d /= tmp;
	
	if(abs(d) < abs(*dist))
	{																					
		vec[0] = d * axb[0]/tmp;
		vec[1] = d * axb[1]/tmp;
		vec[2] = d * axb[2]/tmp;

		*dist = c->direction * d;
	}
	
	//calculate distance to face with normal -axc

	d = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*axc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*axc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*axc[2];
	d *= -1.;
	tmp = sqrt( axc[0]*axc[0] + axc[1]*axc[1] + axc[2]*axc[2] );
	d /= tmp;

	if(abs(d) < abs(*dist))
	{																						
		vec[0] = -d * axc[0]/tmp;
		vec[1] = -d * axc[1]/tmp;
		vec[2] = -d * axc[2]/tmp;

		*dist = c->direction * d;
	}
																		
	//calculate distance to face with normal bxc

	d = (ppos[0]-c->pos[0]-c->a[0]-c->b[0]-c->c[0])*bxc[0] + (ppos[1]-c->pos[1]-c->a[1]-c->b[1]-c->c[1])*bxc[1] + (ppos[2]-c->pos[2]-c->a[2]-c->b[2]-c->c[2])*bxc[2];
	tmp = sqrt( bxc[0]*bxc[0] + bxc[1]*bxc[1] + bxc[2]*bxc[2] );
	d /= tmp;

	if(abs(d) < abs(*dist))
	{																						
		vec[0] = d * bxc[0]/tmp;
		vec[1] = d * bxc[1]/tmp;
		vec[2] = d * bxc[2]/tmp;
	
		*dist = c->direction * d;
	}
}

void calculate_pore_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_pore *c, double *dist, double *vec)
{
  int i;
  double c_dist[3];           /* cartesian distance from pore center */
  double z , r;             /* cylindrical coordinates, coordinate system parallel to
                               pore, origin at pore centera */
  double z_vec[3], r_vec[3]; /* cartesian vectors that correspond to these coordinates */
  double e_z[3], e_r[3];    /* unit vectors in the cylindrical coordinate system */
  /* helper variables, for performance reasons should the be move the the
   * constraint struct*/
     double slope, z_left, z_right;
  /* and another helper that is hopefully optmized out */
     double norm;
     double c1_r, c1_z, c2_r, c2_z;
     double cone_vector_r, cone_vector_z, p1_r, p1_z, dist_vector_z, dist_vector_r, temp;

     
     slope = (c->rad_right - c->rad_left)/2./(c->length-c->smoothing_radius);

  /* compute the position relative to the center of the pore */
  for(i=0;i<3;i++) {
    c_dist[i] = ppos[i] - c->pos[i]; 
  } 
  
  /* compute the component parallel to the pore axis */
  z =0.; 
  for(i=0;i<3;i++) {
    z += (c_dist[i] * c->axis[i]);
  }
  
  /* decompose the position into parallel and perpendicular to the axis */
  r = 0.;
  for(i=0;i<3;i++) {
    z_vec[i] = z * c->axis[i]; 
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
    e_z[i] = c->axis[i];
  norm = 0;
  for(i=0;i<3;i++) 
    norm += r_vec[i]*r_vec[i]; 
  norm = sqrt(norm);
  for(i=0;i<3;i++) 
    e_r[i] = r_vec[i]/norm;
  
  /* c?_r/z and are the centers of the circles that are used to smooth 
   * the entrance of the pore in cylindrical coordinates*/
  c1_z = - (c->length - c->smoothing_radius);
  c2_z = + (c->length - c->smoothing_radius);
  z_left = c1_z - sign(slope) * sqrt(slope*slope/(1+slope*slope))*c->smoothing_radius;
  z_right = c2_z + sign(slope) * sqrt(slope*slope/(1+slope*slope))*c->smoothing_radius;

  c1_r = c->rad_left + slope * ( z_left + c->length ) +
      sqrt( c->smoothing_radius * c->smoothing_radius  - SQR( z_left - c1_z ) );
  c2_r = c->rad_left + slope * ( z_right + c->length ) +
      sqrt( c->smoothing_radius * c->smoothing_radius  - SQR( z_right - c2_z ) );
  c1_r = c->rad_left+c->smoothing_radius;
  c2_r = c->rad_right+c->smoothing_radius;
 
  /* Check if we are in the region of the left wall */
  if (( (r >= c1_r) && (z <= c1_z) ) || ( ( z <= 0 ) && (r>=max(c1_r, c2_r)))) {
    dist_vector_z=-z - c->length;
    dist_vector_r=0;
    *dist = -z - c->length;
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return;
  }
  /* Check if we are in the region of the right wall */
  if (( (r >= c2_r) && (z <= c2_z) ) || ( ( z >= 0 ) && (r>=max(c1_r, c2_r)))) {
    dist_vector_z=-z + c->length;
    dist_vector_r=0;
    *dist = +z - c->length;
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return;
  }

  /* check if the particle should feel the smoothed ends or the middle of the pore */
  /* calculate aufpunkt in z direction first.   */

  /* the distance of the particle from the pore cylinder/cone calculated by projection on the
   * cone normal. Should be > 0 if particle is inside the pore */

  cone_vector_z=1/sqrt(1+slope*slope);
  cone_vector_r=slope/sqrt(1+slope*slope);

  p1_r = c1_r+ ( (r-c1_r)*cone_vector_r + (z-c1_z)*cone_vector_z) * cone_vector_r;
  p1_z = c1_z+ ( (r-c1_r)*cone_vector_r + (z-c1_z)*cone_vector_z) * cone_vector_z;

  dist_vector_r = p1_r-r;
  dist_vector_z = p1_z-z;

  if ( p1_z>=c1_z && p1_z<=c2_z ) {
    if ( dist_vector_r <= 0  ) {
      if (z<0) {
        dist_vector_z=-z - c->length;
        dist_vector_r=0;
        *dist = -z - c->length;
        for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
        return;
      } else {
        dist_vector_z=-z + c->length;
        dist_vector_r=0;
        *dist = +z - c->length;
        for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
        return;
      }
    }
    temp=sqrt( dist_vector_r*dist_vector_r + dist_vector_z*dist_vector_z );
    *dist=temp-c->smoothing_radius;
    dist_vector_r-=dist_vector_r/temp*c->smoothing_radius;
    dist_vector_z-=dist_vector_z/temp*c->smoothing_radius;
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return;
  }


  /* Check if we are in the range of the left smoothing circle */
  if (p1_z <= c1_z ) {
    /* distance from the smoothing center */
    norm = sqrt( (z - c1_z)*(z - c1_z) + (r - c1_r)*(r - c1_r) );
    *dist = norm - c->smoothing_radius;
    dist_vector_r=(c->smoothing_radius/norm -1)*(r - c1_r);
    dist_vector_z=(c->smoothing_radius/norm - 1)*(z - c1_z);
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return;
  }
  /* Check if we are in the range of the right smoothing circle */
  if (p1_z >= c2_z ) {
    norm = sqrt( (z - c2_z)*(z - c2_z) + (r - c2_r)*(r - c2_r) );
    *dist = norm - c->smoothing_radius;
    dist_vector_r=(c->smoothing_radius/norm -1)*(r - c2_r);
    dist_vector_z=(c->smoothing_radius/norm - 1)*(z - c2_z);
    for (i=0; i<3; i++) vec[i]=-dist_vector_r*e_r[i] - dist_vector_z*e_z[i];
    return;
  }
  exit(printf("should never be reached, z %f, r%f\n",z, r));
}

void calculate_plane_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_plane *c, double *dist, double *vec)
{
  int i;
  double c_dist_sqr,c_dist;
  
  c_dist_sqr=0.0;
  for(i=0;i<3;i++) {
    if(c->pos[i] >= 0) {
      vec[i] = c->pos[i] - ppos[i];
      c_dist_sqr += SQR(vec[i]);
    }else{
      vec[i] = 0.0;
      c_dist += SQR(vec[i]);
    }
  }
  c_dist = sqrt(c_dist_sqr);
  *dist = c_dist;

  
  for(i=0;i<3;i++) {
    vec[i] *= -1;
  }
}

void add_rod_force(Particle *p1, double ppos[3], Particle *c_p, Constraint_rod *c)
{
#ifdef ELECTROSTATICS
  int i;
  double fac, vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = ppos[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->lambda != 0.0) {
    fac = 2*coulomb.prefactor*c->lambda*p1->p.q/c_dist_2;
    p1->f.f[0]  += fac*vec[0];
    p1->f.f[1]  += fac*vec[1];
    c_p->f.f[0] -= fac*vec[0];
    c_p->f.f[1] -= fac*vec[1];
  }
#endif
}

double rod_energy(Particle *p1, double ppos[3], Particle *c_p, Constraint_rod *c)
{
#ifdef ELECTROSTATICS
  int i;
  double vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = ppos[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->lambda != 0.0)
    return coulomb.prefactor*p1->p.q*c->lambda*(-log(c_dist_2*SQR(box_l_i[2])) + 2*(M_LN2 - C_GAMMA));
#endif
  return 0;
}

void add_plate_force(Particle *p1, double ppos[3], Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  double f;

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0) {
    f = 2*M_PI*coulomb.prefactor*c->sigma*p1->p.q;
    if (ppos[2] < c->pos)
      f = -f;
    p1->f.f[2]  += f;
    c_p->f.f[2] -= f;
  }
#endif
}

double plate_energy(Particle *p1, double ppos[3], Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0)
    return -2*M_PI*coulomb.prefactor*c->sigma*p1->p.q*fabs(ppos[2] - c->pos);
#endif
  return 0;
}

//ER
void add_ext_magn_field_force(Particle *p1, Constraint_ext_magn_field *c)
{
#ifdef ROTATION
#ifdef DIPOLES
  p1->f.torque[0] += p1->r.dip[1]*c->ext_magn_field[2]-p1->r.dip[2]*c->ext_magn_field[1];
  p1->f.torque[1] += p1->r.dip[2]*c->ext_magn_field[0]-p1->r.dip[0]*c->ext_magn_field[2];
  p1->f.torque[2] += p1->r.dip[0]*c->ext_magn_field[1]-p1->r.dip[1]*c->ext_magn_field[0];
#endif
#endif
}

double ext_magn_field_energy(Particle *p1, Constraint_ext_magn_field *c)
{
#ifdef DIPOLES
     return -1.0 * scalar(c->ext_magn_field,p1->r.dip);
#endif
  return 0;
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

void add_constraints_forces(Particle *p1)
{
  if (n_constraints==0)
   return;
  int n, j;
  double dist, vec[3], force[3], torque1[3], torque2[3];

  IA_parameters *ia_params;
  char *errtxt;
  double folded_pos[3];
  int img[3];

  /* fold the coordinate[2] of the particle */
  memcpy(folded_pos, p1->r.p, 3*sizeof(double));
  memcpy(img, p1->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);

  for(n=0;n<n_constraints;n++) {
    ia_params=get_ia_param(p1->p.type, (&constraints[n].part_rep)->p.type);
    dist=0.;
    for (j = 0; j < 3; j++) {
      force[j] = 0;
#ifdef ROTATION
      torque1[j] = torque2[j] = 0;
#endif
    }

    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      if(checkIfInteraction(ia_params)) {
	calculate_wall_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.wal.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.wal.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.wal.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{061 wall constraint %d violated by particle %d} ", n, p1->p.identity);
    }
	}
      }
      break;

    case CONSTRAINT_SPH:
      if(checkIfInteraction(ia_params)) {
	calculate_sphere_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.sph.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.sph.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.sph.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{062 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
    }
	}
      }
      break;
    
    case CONSTRAINT_CYL: 
      if(checkIfInteraction(ia_params)) {
	calculate_cylinder_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.cyl, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.cyl.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.cyl.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.cyl.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
    }
        }
      }
      break;
    
    case CONSTRAINT_RHOMBOID: 
      if(checkIfInteraction(ia_params)) {
	calculate_rhomboid_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.rhomboid, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.rhomboid.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.rhomboid.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.rhomboid.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 rhomboid constraint %d violated by particle %d} ", n, p1->p.identity);
    }
        }
      }
      break;
	
    case CONSTRAINT_MAZE: 
      if(checkIfInteraction(ia_params)) {
	calculate_maze_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.maze.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{064 maze constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_PORE: 
      if(checkIfInteraction(ia_params)) {
	calculate_pore_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist, vec); 
	if ( dist >= 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else {
    if(constraints[n].c.pore.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.pore.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 pore constraint %d violated by particle %d} ", n, p1->p.identity);
        }
      }
      }
      break;
	
      /* electrostatic "constraints" */
    case CONSTRAINT_ROD:
      add_rod_force(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.rod);
      break;

    case CONSTRAINT_PLATE:
      add_plate_force(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.plate);
      break;
      
#ifdef DIPOLES
    case CONSTRAINT_EXT_MAGN_FIELD:
      add_ext_magn_field_force(p1, &constraints[n].c.emfield);
      break;
#endif
    
    case CONSTRAINT_PLANE:
     if(checkIfInteraction(ia_params)) {
	calculate_plane_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.plane, &dist, vec); 
	if (dist > 0) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				       ia_params,vec,dist,dist*dist, force,
				       torque1, torque2);
#ifdef TUNABLE_SLIP
	    add_tunable_slip_pair_force(p1, &constraints[n].part_rep,ia_params,vec,dist,force);
#endif
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 plane constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
    }
    for (j = 0; j < 3; j++) {
      p1->f.f[j] += force[j];
      constraints[n].part_rep.f.f[j] -= force[j];
#ifdef ROTATION
      p1->f.torque[j] += torque1[j];
      constraints[n].part_rep.f.torque[j] += torque2[j];
#endif
    }
  }
}

double add_constraints_energy(Particle *p1)
{
  int n, type;
  double dist, vec[3];
  double nonbonded_en, coulomb_en, magnetic_en;
  IA_parameters *ia_params;
  char *errtxt;
  double folded_pos[3];
  int img[3];

  /* fold the coordinate[2] of the particle */
  memcpy(folded_pos, p1->r.p, 3*sizeof(double));
  memcpy(img, p1->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
  for(n=0;n<n_constraints;n++) { 
    ia_params = get_ia_param(p1->p.type, (&constraints[n].part_rep)->p.type);
    nonbonded_en = 0.;
    coulomb_en   = 0.;
    magnetic_en = 0.;

    dist=0.;
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      if(checkIfInteraction(ia_params)) {
	calculate_wall_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else if ( dist <= 0 && constraints[n].c.wal.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{065 wall constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_SPH: 
      if(checkIfInteraction(ia_params)) {
	calculate_sphere_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else if ( dist <= 0 && constraints[n].c.sph.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{066 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_CYL: 
      if(checkIfInteraction(ia_params)) {
	calculate_cylinder_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.cyl, &dist , vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);

	}
	else if ( dist <= 0 && constraints[n].c.cyl.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_RHOMBOID: 
      if(checkIfInteraction(ia_params)) {
	calculate_rhomboid_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.rhomboid, &dist , vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);

	}
	else if ( dist <= 0 && constraints[n].c.rhomboid.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_MAZE: 
      if(checkIfInteraction(ia_params)) {
	calculate_maze_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else if ( dist <= 0 && constraints[n].c.maze.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{068 maze constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_PORE: 
      if(checkIfInteraction(ia_params)) {
	calculate_pore_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist , vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);

	}
	else {
	  errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 pore constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_ROD:
      coulomb_en = rod_energy(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.rod);
      break;

    case CONSTRAINT_PLATE:
      coulomb_en = plate_energy(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.plate);
      break;
    case CONSTRAINT_EXT_MAGN_FIELD:
      magnetic_en = ext_magn_field_energy(p1, &constraints[n].c.emfield);
      break;
    }

    if (energy.n_coulomb > 0)
      energy.coulomb[0] += coulomb_en;
    
    if (energy.n_dipolar > 0)
      energy.dipolar[0] += magnetic_en;

    type = (&constraints[n].part_rep)->p.type;
    if (type >= 0)
      *obsstat_nonbonded(&energy, p1->p.type, type) += nonbonded_en;
  }
  return 0.;
}

#endif

