// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef NEMD_H
#define NEMD_H
/** \file nemd.h

    This file contains the implementation of the NEMD (Non Equilibrium
    Molecular Dynamics) algorithm. It allows one to shear a system
    with help of an unphysical momentum change in two slabs in the
    system.

    Usage:
    <ul>
    <li> \verbatim nemd \endverbatim        Return nemd status.
    <li> \verbatim nemd off \endverbatim    Turn off nemd.
    <li> \verbatim nemd exchange <INT n_slabs> <INT n_exchange> \endverbatim
         enable method exchange momentum.
    <li> \verbatim nemd shearrate <INT n_slabs> <DOUBLE shearrate> \endverbatim
         enable method with fixed shear rate.
    <li> \verbatim nemd profile \endverbatim    Return the velocity profile.
    <li> \verbatim nemd viscosity \endverbatim  Return the viscosity.
    </ul>

    Notes:
    <ul>
    <li> The slabs are oriented perpendicualr to the z direction.
    <li> The velocity profile is generated in x direction.
    <li> Do not use other special features like part fix or 
         constraints inside the top and middle slab.
    <li> Use only with a DPD thermostat or in an NVE ensemble.
    </ul>
    
    Methods:

    Both methods devide the simulation box into n_slab slabs (x-y
    planes) in z direction and apply a shear in x direction. The shear
    is aplied on the top and the middle layer.
    
    <ul>
    <li> <b>Exchange</b>: At each step find the n_exchange particles
    with the largest velocity x-component (top slab negative velocity
    component, middle slab positiv velocity component) in the top and
    middle slab. Exchange the x component of their velocities.
    <li> <b>Shearrate</b>: Calculate the mean velocity of the top and
    middle slab that produces the desired shear rate: mean_velocity =
    (+/-) shear rate*box_l[2]/4.0. During integration calculate the
    actual mean velocities of these slabs and add the difference
    between the actual mean velocity and the desired mean velocity to
    the velocity of the particles in the slab.
    </ul>

    Results:

    During the simulation the velocity profile along the slabs as well
    as the viscosity is calculated.
    The viscosity is calculated via:
    \f[ \eta = \frac{F}{\dot{\gamma} L_x L_y} \f]
    F is the mean force (momentum transfer per unit time) acting on
    the slab and \f[ L_x L_y \f] is the area of the slab.

    <b>Responsible:</b>
    <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 */

#include <tcl.h>
#include <mpi_wrap.h>
#include <stdio.h>
#include "particle_data.h"
#include "grid.h"
#include "global.h"
#include "parser.h"


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

/** tcl procedure for nemd steering.
    USAGE: nemd \<n_slabs\> \<n_exchange\>   
    see also \ref tcl_nemd
*/
int nemd(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

#ifdef NEMD 

/** Change momentum in top an middle slab. How this is done depends on
    nemd method. */
void nemd_change_momentum();

/** Store the mean value of the velocity x-component of all slabs in
    \ref Nemd::velocity_profile. */ 
void nemd_store_velocity_profile();

/** Store the x-componenet of the velocity of particle part into the
    nemd data structure \ref nemddata. */
MDINLINE void nemd_get_velocity(Particle part) 
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

MDINLINE void nemd_add_velocity (Particle *part) {
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
