// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef NEMD_H
#define NEMD_H
/** \file nemd.h

    This file contains the implementation of the NEMD (Non Equilibrium
    Molecular Dynamics) algorithm. It allows one to shear a system
    with help of an unphysical momentum exchange between two slabs in
    the system.

    The slabs are oriented perpendicualr to the z direction!

    <b>Responsible:</b>
    <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 */

#include <tcl.h>
#include <mpi.h>
#include <stdio.h>
#include "particle_data.h"
#include "grid.h"
#include "global.h"
#include "parser.h"
#include "integrate.h"

#ifdef NEMD

#define NEMD_METHOD_NOTSET    0
#define NEMD_METHOD_EXCHANGE  1
#define NEMD_METHOD_SHEARRATE 2


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

  /** used method to create the shear */
  int method;

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

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure for nemd steering.
    USAGE: nemd \<n_slabs\> \<n_exchange\>   
    see also \ref tcl_nemd
*/
int nemd(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);
/** Initialize nemd data structure \ref nemddata for the values
    n_slabs and n_exchange. If n_slabs is 0 it frees the nemd data
    structure \ref nemddata and turns of the nemd part of the
    integrator. Called by \ref nemd. */
void nemd_init(int n_slabs, int n_exchange, double shear_rate);

/** Free an old nemd data structure \ref nemddata and set its value to
    the initialization values. Called by nemd_init. */
void nemd_free();

/** Exchange momentum between top an middle slab. */
void nemd_change_momentum();

/** Store the mean value of the velocity x-component of all slabs in
    \ref nemddata::velocity_profile. */ 
void nemd_store_velocity_profile();

/** Print out velocity profile to tcl */
int nemd_print_profile(Tcl_Interp *interp);

/** Store the x-componenet of the velocity of particle part into the
    nemd data structure \ref nemddata. */
MDINLINE void nemd_store_velocity(Particle part) {
#ifdef NEMD 
  int slab_num, i;
  Slab *slab;

  if(nemddata.n_slabs == -1 ) return;
  if(nemddata.method == NEMD_METHOD_EXCHANGE && nemddata.n_exchange==0 ) return;
  //  INTEG_TRACE(fprintf(stderr,"%d: nemd_store_velocity: Part %d in %d slabs\n",this_node,part.p.identity,nemddata.n_slabs));


  slab_num = part.r.p[2]*nemddata.invthickness;
  /* THIS IS A SINGLE NODE CORRECTION !!! */
  if(slab_num == nemddata.n_slabs) slab_num = 0;
  if(part.r.p[2] < 0.0) slab_num = nemddata.n_slabs-1;

  slab = &nemddata.slab[slab_num];
  
  if(nemddata.method == NEMD_METHOD_EXCHANGE) {
    /* look for largest negative x componenet */
    if(slab_num == nemddata.top_slab) {
    
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
      /* look for largest positive x componenet */
    } else if(slab_num == nemddata.mid_slab) {
      
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
    /* Add the velocities to mean velocity of that slab */
    slab->v_mean          += part.m.v[0];
    slab->n_parts_in_slab ++;
  }
  else if(nemddata.method == NEMD_METHOD_SHEARRATE) {
    if(slab_num == nemddata.top_slab || slab_num == nemddata.mid_slab) {
      part.m.v[0] += slab->vel_diff;
    } 
    /* Add the velocities = velocity+force to mean velocity of that slab */
    slab->v_mean          += part.m.v[0] + part.f.f[0];
    slab->n_parts_in_slab ++;

  }
#endif
}


MDINLINE void nemd_add_velocity (Particle *part) {
#ifdef NEMD 
  int slab_num;
  Slab *slab;

  if( !(nemddata.method == NEMD_METHOD_SHEARRATE) ) return;

  slab_num = part->r.p[2]*nemddata.invthickness;
  slab = &nemddata.slab[slab_num];

  if(slab_num == nemddata.top_slab || slab_num == nemddata.mid_slab) {
    part->f.f[0] += slab->vel_diff;
  }
#endif
}


/*@}*/


#endif
