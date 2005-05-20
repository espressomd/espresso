// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#ifndef CONSTRAINT_H
#define CONSTRAINT_H
#include "statistics.h"
#include "energy.h"
#include "grid.h"
#include "errorhandling.h"

/** \file constraint.h
 *  Routines for handling of constraints. Implemented are walls, cylinders and spheres.
 *  Only activ if you define CONSTRAINTS in \ref config.h.
 *  see also \ref interaction_data.h
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
 */

#ifdef CONSTRAINTS
#ifndef LENNARD_JONES
#error **********************************
#error
#error CONSTRAINTS requires LENNARD_JONES
#error
#error **********************************
#else

#define C_GAMMA   0.57721566490153286060651209008

MDINLINE void calculate_wall_dist(Particle *p1, Particle *c_p, Constraint_wall *c, double *dist, double *vec)
{
  int i;
  IA_parameters *ia_params;

  ia_params=get_ia_param(p1->p.type, c_p->p.type);

  *dist = -c->d;
  for(i=0;i<3;i++) *dist += p1->r.p[i]*c->n[i];
  
  for(i=0;i<3;i++) vec[i] = c->n[i] * *dist;
  
}


MDINLINE void calculate_sphere_dist(Particle *p1, Particle *c_p, Constraint_sphere *c, double *dist, double *vec)
{
  int i;
  double fac,  c_dist;
  IA_parameters *ia_params;
 
  ia_params=get_ia_param(p1->p.type, c_p->p.type);
 
  c_dist=0.0;
  for(i=0;i<3;i++) {
    vec[i] = c->pos[i] - p1->r.p[i];
    c_dist += SQR(vec[i]);
  }
  c_dist = sqrt(c_dist);
  
  if ( c->direction == -1 ) {
  /* apply force towards inside the sphere */
    *dist = c->rad - c_dist;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= fac;}
  else  {
   /* apply force towards outside the sphere */
    *dist = c_dist - c->rad;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= -fac;
  }
}


MDINLINE void calculate_maze_dist(Particle *p1, Particle *c_p, Constraint_maze *c, double *dist, double *vec)
{
  int i,min_axis,cursph[3],dim;
  float diasph,fac,c_dist,sph_dist,cyl_dist,temp_dis;
  float sph_vec[3],cyl_vec[3];
  IA_parameters *ia_params;
 
  ia_params = get_ia_param(p1->p.type, c_p->p.type);
  dim=(int) c->dim;
  diasph = box_l[0]/c->nsphere;
  
  /* First determine the distance to the sphere */
  c_dist=0.0;
  for(i=0;i<3;i++) {
    cursph[i] = (int) (p1->r.p[i]/diasph);
    sph_vec[i] = (cursph[i]+0.5) * diasph  - (p1->r.p[i]);
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

MDINLINE void calculate_cylinder_dist(Particle *p1, Particle *c_p, Constraint_cylinder *c, double *dist, double *vec)
{
  int i;
  double d_per,d_par,d_real,d_per_vec[3],d_par_vec[3],d_real_vec[3];
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->p.type, c_p->p.type);


  d_real = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = p1->r.p[i] - c->pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);
    
  d_par=0.;
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * c->axis[i]);
  }
    
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * c->axis[i] ;
    d_per_vec[i] = p1->r.p[i] - (c->pos[i] + d_par_vec[i]) ;
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

MDINLINE void add_rod_force(Particle *p1, Particle *c_p, Constraint_rod *c)
{
#ifdef ELECTROSTATICS
  int i;
  double fac, vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = p1->r.p[i] - c->pos[i];
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

MDINLINE double rod_energy(Particle *p1, Particle *c_p, Constraint_rod *c)
{
#ifdef ELECTROSTATICS
  int i;
  double vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = p1->r.p[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->lambda != 0.0)
    return coulomb.prefactor*p1->p.q*c->lambda*(-log(c_dist_2*SQR(box_l_i[2])) + 2*(M_LN2 - C_GAMMA));
#endif
  return 0;
}

MDINLINE void add_plate_force(Particle *p1, Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  double f;
  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0) {
    f = 2*M_PI*coulomb.prefactor*c->sigma*p1->p.q;
    if (p1->r.p[2] < c->pos)
      f = -f;
    p1->f.f[2]  += f;
    c_p->f.f[2] -= f;
  }
#endif
}

MDINLINE double plate_energy(Particle *p1, Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0)
    return 2*M_PI*coulomb.prefactor*c->sigma*p1->p.q*fabs(p1->r.p[2] < c->pos);
#endif
  return 0;
}

MDINLINE void add_constraints_forces(Particle *p1)
{
  int n, j;
  double dist, vec[3], force[3];
  IA_parameters *ia_params;
  char *errtxt;

  for(n=0;n<n_constraints;n++) {
    ia_params=get_ia_param(p1->p.type, (&constraints[n].part_rep)->p.type);
    dist=0.;
    for (j = 0; j < 3; j++)
      force[j] = 0;
	
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      if(ia_params->LJ_cut > 0. ) {
	calculate_wall_dist(p1, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
	if (dist > 0) {
	  add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist, force);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{061 wall constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_SPH:
      if(ia_params->LJ_cut > 0. ) {
	calculate_sphere_dist(p1, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if (dist > 0) {
	  add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist, force);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{062 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
    
    case CONSTRAINT_CYL: 
      if(ia_params->LJ_cut > 0. ) {
	calculate_cylinder_dist(p1, &constraints[n].part_rep, &constraints[n].c.cyl, &dist, vec); 
	if ( dist > 0 ) {
	  add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist, force);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
        }
      }
      break;
	
    case CONSTRAINT_MAZE: 
      calculate_maze_dist(p1, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
      if (dist > 0) {
	add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist, force);
      }
      else {
	errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt, "{064 maze constraint %d violated by particle %d} ", n, p1->p.identity);
      }
      break;

      /* electrostatic "constraints" */
    case CONSTRAINT_ROD:
      add_rod_force(p1, &constraints[n].part_rep, &constraints[n].c.rod);
      break;

    case CONSTRAINT_PLATE:
      add_plate_force(p1, &constraints[n].part_rep, &constraints[n].c.plate);
      break;
    }

    for (j = 0; j < 3; j++) {
      p1->f.f[j] += force[j];
      constraints[n].part_rep.f.f[j] -= force[j];
    }
  }
}

MDINLINE double add_constraints_energy(Particle *p1)
{
  int n, type;
  double dist, vec[3];
  double lj_en, coulomb_en;
  IA_parameters *ia_params;
  char *errtxt;

  for(n=0;n<n_constraints;n++) { 
    ia_params = get_ia_param(p1->p.type, (&constraints[n].part_rep)->p.type);
    lj_en      = 0;
    coulomb_en = 0;

    dist=0.;
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      if(ia_params->LJ_cut > 0. ) {
	calculate_wall_dist(p1, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
	if (dist > 0)
	  lj_en = lj_pair_energy(p1, &constraints[n].part_rep, ia_params, vec, dist);
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{065 wall constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_SPH: 
      if(ia_params->LJ_cut > 0. ) {
	calculate_sphere_dist(p1, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if (dist > 0) {
	  lj_en = lj_pair_energy(p1, &constraints[n].part_rep, ia_params, vec, dist);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{066 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_CYL: 
      if(ia_params->LJ_cut > 0. ) {
	calculate_cylinder_dist(p1, &constraints[n].part_rep, &constraints[n].c.cyl, &dist , vec); 
	if ( dist > 0 ) {
	  lj_en = lj_pair_energy(&constraints[n].part_rep, p1, ia_params, vec, dist);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_MAZE: 
      calculate_maze_dist(p1, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
      if (dist > 0) {
	lj_en = lj_pair_energy(&constraints[n].part_rep, p1, ia_params, vec, dist);
      }
      else {
	errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt, "{068 maze constraint %d violated by particle %d} ", n, p1->p.identity);
      }
      break;

    case CONSTRAINT_ROD:
      coulomb_en = rod_energy(p1, &constraints[n].part_rep, &constraints[n].c.rod);
      break;

    case CONSTRAINT_PLATE:
      coulomb_en = plate_energy(p1, &constraints[n].part_rep, &constraints[n].c.plate);
      break;
    }

    if (energy.n_coulomb > 0)
      energy.coulomb[0] += coulomb_en;

    type = (&constraints[n].part_rep)->p.type;
    if (type >= 0)
      *obsstat_nonbonded(&energy, p1->p.type, type) += lj_en;
  }
  return 0.;
}

MDINLINE void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < n_constraints; n++)
    for (i = 0; i < 3; i++)
      constraints[n].part_rep.f.f[i] = 0;
}
#endif

#endif
#endif
