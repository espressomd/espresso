/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#ifndef NEMD_H
#define NEMD_H
/** \file nemd.hpp

    This file contains the implementation of the NEMD (Non Equilibrium
    Molecular Dynamics) algorithm. It allows one to shear a system
    with help of an unphysical momentum change in two slabs in the
    system.
 */

#include "particle_data.hpp"

#define NEMD_METHOD_OFF       0
#define NEMD_METHOD_EXCHANGE  1
#define NEMD_METHOD_SHEARRATE 2

#ifdef NEMD

/************************************************/
/** \name Data Types */
/************************************************/
/*@{*/
/** Data structure describing a slab and the velocities occuring their in. */
typedef struct {

  /* velocity profile */
  /** mean velocity in slab */
  double v_mean;
  /** Number of particles in slab */
  int n_parts_in_slab;

  /* the n_exchange fastest particles */
  /** identity list of n_exchange fastest particles */
  int *fastest;
  /** number of particles allready stored */
  int n_fastest;
  /** minimal velocity stored in fastest*/
  double v_min;
  /** index of that particle in fastest */
  int ind_min;

  /* difference to desired mean velocity in slab */
  double vel_diff;
} Slab;

/** Structure containing the NEMD informations */
typedef struct {
  /** number of slabs */
  int n_slabs;
  /** index of middle slab */
  int mid_slab;
  /** index of top slab */
  int top_slab;

  /** thickness of a slab */
  double thickness;
  /** invers thickness of a slab */
  double invthickness;

  /** Velocity profile */
  double *velocity_profile;
  /** Number of profiles stored */
  int profile_norm;

 /* method 1: exchange momentum */
  /** Number of particle momentums to exchange */
  int n_exchange;
  /** Array containing the slabs */
  Slab *slab;

  /* method 2: shear rate */
  /** Desired shear rate */
  double shear_rate;
  /** Corresponding velocity of top and mid slab */
  double slab_vel;

  /* momentum added/exchanged */
  double momentum;
  /* numbe of momentum exchanges */
  int momentum_norm;
} Nemd;
/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Structure containing the nemd relevant information */
extern Nemd nemddata;
/*@}*/
#endif

extern int nemd_method;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

#ifdef NEMD 

int nemd_free(void);

void nemd_init(int n_slabs, int n_exchange, double shear_rate);


/** Change momentum in top an middle slab. How this is done depends on
    nemd method. */
void nemd_change_momentum();

/** Store the mean value of the velocity x-component of all slabs in
    \ref Nemd::velocity_profile. */ 
void nemd_store_velocity_profile();

/** Store the x-componenet of the velocity of particle part into the
    nemd data structure \ref nemddata. */
inline void nemd_get_velocity(Particle part) 
{
  int  slab_num, i;
  Slab *slab;

  if(nemd_method == NEMD_METHOD_OFF ) return;
  if(nemd_method == NEMD_METHOD_EXCHANGE && nemddata.n_exchange==0 ) return;

  /* calculate slab_num */
  slab_num = part.r.p[2]*nemddata.invthickness;
  /* THIS IS A SINGLE NODE CORRECTION !!! */
  if(slab_num == nemddata.n_slabs) slab_num = 0;
  if(part.r.p[2] < 0.0) slab_num = nemddata.n_slabs-1;
  slab = &nemddata.slab[slab_num];
  
  /* Add velocity to mean velocity of that slab */
  slab->v_mean          += part.m.v[0];
  slab->n_parts_in_slab ++;

  /* Collect the n_exchange particles with largest velocities for top
     and middle slab if nemd method is exchange */
  if(nemd_method == NEMD_METHOD_EXCHANGE) {
    if(slab_num == nemddata.top_slab) {
      /* top slab: look for largest negative x componenet */
   
      /* first fill fastest array as particles com in */
      if(slab->n_fastest < nemddata.n_exchange) {     
	slab->fastest[slab->n_fastest] = part.p.identity;
	if( part.m.v[0] > slab->v_min ) {
	  slab->v_min   = part.m.v[0];
	  slab->ind_min = slab->n_fastest; 
	}
	slab->n_fastest++;
      } 
      /* then replace the slowest with a new one if necessary... */  
      else if ( part.m.v[0] < slab->v_min ) {
	slab->fastest[slab->ind_min] = part.p.identity;
	/* ... and find again the slowest one now */
	slab->v_min = local_particles[slab->fastest[0]]->m.v[0];	  
	slab->ind_min = 0;
	for(i=1;i<slab->n_fastest;i++) 
	  if(local_particles[slab->fastest[i]]->m.v[0] > slab->v_min) {
	    slab->v_min   = local_particles[slab->fastest[i]]->m.v[0];
	    slab->ind_min = i;
	  }       
      } 
    } else if(slab_num == nemddata.mid_slab) {
      /* midlle slab: look for largest positive x componenet */

      /* first fill fastest array as particles com in */
      if(slab->n_fastest < nemddata.n_exchange) {
	slab->fastest[slab->n_fastest] = part.p.identity;
	if( part.m.v[0] < slab->v_min ) {
	  slab->v_min   = part.m.v[0];
	  slab->ind_min = slab->n_fastest; 
	}
	slab->n_fastest++;
      } 
      /* then replace the slowest with a new one if necessary... */  
      else if ( part.m.v[0] > slab->v_min) {
	slab->fastest[slab->ind_min] = part.p.identity;
	/* ... and find again the slowest one now */
	slab->v_min = local_particles[slab->fastest[0]]->m.v[0];
	slab->ind_min = 0;
	for(i=1;i<slab->n_fastest;i++) 
	  if(local_particles[slab->fastest[i]]->m.v[0] < slab->v_min) {
	    slab->v_min   = local_particles[slab->fastest[i]]->m.v[0];
	    slab->ind_min = i;
	  }
      }    
    } 
  }
}

inline void nemd_add_velocity (Particle *part) {
  if( !(nemd_method == NEMD_METHOD_SHEARRATE) ) {
    return;
  } else {
    int slab_num;
    Slab *slab;

    slab_num = part->r.p[2]*nemddata.invthickness;
    slab = &nemddata.slab[slab_num];
    
    if(slab_num == nemddata.top_slab || slab_num == nemddata.mid_slab) {
      part->m.v[0] += slab->vel_diff;
    }
  }
}

/* endif NEMD */
#endif

/*@}*/

/* endif NEMD_H */
#endif
