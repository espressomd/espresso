/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef CONSTRAINT_H
#define CONSTRAINT_H
#include "statistics.h"
#include "energy.h"
#include "forces.h"
#include "grid.h"
#include "errorhandling.h"
#include "tunable_slip.h"

/** \file constraint.h
 *  Routines for handling of constraints.
 *  Only active if the feature CONSTRAINTS is activated.
 *  see also \ref interaction_data.h
 */

#ifdef CONSTRAINTS

// for the charged rod "constraint"
#define C_GAMMA   0.57721566490153286060651209008

MDINLINE void calculate_wall_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_wall *c, double *dist, double *vec)
{
  int i;

  *dist = -c->d;
  for(i=0;i<3;i++) *dist += ppos[i]*c->n[i];
  
  for(i=0;i<3;i++) vec[i] = c->n[i] * *dist;
  
}


MDINLINE void calculate_sphere_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_sphere *c, double *dist, double *vec)
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


MDINLINE void calculate_maze_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_maze *c, double *dist, double *vec)
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

MDINLINE void calculate_cylinder_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_cylinder *c, double *dist, double *vec)
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

MDINLINE void calculate_pore_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_pore *c, double *dist, double *vec)
{
  int i;
  double d_per, d_per_sh, d_par, d_par_sh, d_real, d_real_2,
    d_per_vec[3], d_par_vec[3], d_real_vec[3];

  /* compute the position relative to the center of the pore */
  d_real_2 = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = ppos[i] - c->pos[i]; 
    d_real_2 += SQR(d_real_vec[i]);  
  } 
  d_real = sqrt(d_real_2);
  
  /* compute the component parallel to the pore axis */
  d_par=0.; 
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * c->axis[i]);
  }
  
  /* decompose the position into parallel and perpendicular to the axis */
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * c->axis[i]; 
    d_per_vec[i] = d_real_vec[i] - d_par_vec[i];
  } 
  /* compute the norm of the above vectors */
  d_per = sqrt(fabs(d_real_2 - SQR(d_par)));
  d_par = fabs(d_par);

  /* compute the distance from the pore (cylinder and wall) surfaces */
  d_par_sh = d_par - c->length;
  d_per_sh = c->rad - d_per;

  if(d_par_sh <= 0.) { /* within the channel OR the constraint violated */
    if(d_per_sh <= 0.) { /* constraint violated */
	    *dist=-1.0;
	    vec[0]=vec[1]=vec[2]=0.0;
    } 
    else {           /* within the channel */
	    *dist = d_per_sh;   
            for (i=0; i<3;i++) 
               vec[i]= -d_per_vec[i] * (d_per_sh /  d_per);
    }
  }
  else { /* out of the channel (aligned or not) ( d_par_sh > 0 ) */
	  if (d_per_sh < 0) /* not aligned, feel the wall */
	  {
             *dist = d_par_sh;
             for (i=0; i<3;i++) 
                 vec[i]= d_par_vec[i] * (d_par_sh /  d_par) ;
	  }
	  else /* aligned, feel the "ring" */
	  {
	     *dist = sqrt(SQR(d_par_sh)+SQR(d_per_sh));
             for (i=0; i<3;i++) 
               vec[i]= -d_per_vec[i] * (d_per_sh /  d_per)  + 
	                d_par_vec[i] * (d_par_sh /  d_par);
	  }
  }
}

MDINLINE void calculate_plane_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_plane *c, double *dist, double *vec)
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

MDINLINE void add_rod_force(Particle *p1, double ppos[3], Particle *c_p, Constraint_rod *c)
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

MDINLINE double rod_energy(Particle *p1, double ppos[3], Particle *c_p, Constraint_rod *c)
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

MDINLINE void add_plate_force(Particle *p1, double ppos[3], Particle *c_p, Constraint_plate *c)
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

MDINLINE double plate_energy(Particle *p1, double ppos[3], Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0)
    return -2*M_PI*coulomb.prefactor*c->sigma*p1->p.q*fabs(ppos[2] - c->pos);
#endif
  return 0;
}

//ER
MDINLINE void add_ext_magn_field_force(Particle *p1, Constraint_ext_magn_field *c)
{
#ifdef ROTATION
#ifdef DIPOLES
  p1->f.torque[0] += p1->r.dip[1]*c->ext_magn_field[2]-p1->r.dip[2]*c->ext_magn_field[1];
  p1->f.torque[1] += p1->r.dip[2]*c->ext_magn_field[0]-p1->r.dip[0]*c->ext_magn_field[2];
  p1->f.torque[2] += p1->r.dip[0]*c->ext_magn_field[1]-p1->r.dip[1]*c->ext_magn_field[0];
#endif
#endif
}

MDINLINE double ext_magn_field_energy(Particle *p1, Constraint_ext_magn_field *c)
{
#ifdef DIPOLES
     return -1.0 * scalar(c->ext_magn_field,p1->r.dip);
#endif
  return 0;
}
//end ER

MDINLINE void add_constraints_forces(Particle *p1)
{
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
	if (dist > 0) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{061 wall constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_SPH:
      if(checkIfInteraction(ia_params)) {
	calculate_sphere_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if (dist > 0) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{062 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
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
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
        }
      }
      break;
	
    case CONSTRAINT_MAZE: 
      if(checkIfInteraction(ia_params)) {
	calculate_maze_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
	if (dist > 0) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{064 maze constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_PORE: 
      if(checkIfInteraction(ia_params)) {
	calculate_pore_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 pore constraint %d violated by particle %d} ", n, p1->p.identity);
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
      
#ifdef MAGNETOSTATICS 
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
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
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

MDINLINE double add_constraints_energy(Particle *p1)
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
	if (dist > 0)
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{065 wall constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_SPH: 
      if(checkIfInteraction(ia_params)) {
	calculate_sphere_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if (dist > 0) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
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
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_MAZE: 
      if(checkIfInteraction(ia_params)) {
	calculate_maze_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
	if (dist > 0) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
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
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
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

MDINLINE void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < n_constraints; n++)
    for (i = 0; i < 3; i++)
      constraints[n].part_rep.f.f[i] = 0;
}
#endif

#endif
