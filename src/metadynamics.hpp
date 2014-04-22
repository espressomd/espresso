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

#ifndef METADYNAMICS_H
#define METADYNAMICS_H

#include <cstring>
#include <cmath>
#include "particle_data.hpp"
#include "utils.hpp"

/** \file metadynamics.hpp 
*
*  This file contains routines to perform metadynamics.  Right now, the
*  reaction coordinate is defined between two particles. Note that these
*  particles can be virtual sites, in order to handle molecules.
*
*  - set metadynamics options 
*  - initialize bias forces and free energy profiles 
*  - calculate reaction coordinate for each integration step 
*  - apply bias force on particles
*/

#ifdef METADYNAMICS

/** \name Metadynamics switches */
/********************************/
/* off */
#define META_OFF      0
/* distance between two particles */
#define META_DIST     1
/* relative height (z coord) of meta_pid1 with respect to meta_pid2 
* Example: measure height of particle (pid1) with respect to interface (pid2).
*/
#define META_REL_Z   2

/**********************************
* exported variables
**********************************/

/** Flag - turn metadynamics on */
extern int    meta_switch;
/** pid of particle 1 */
extern int    meta_pid1;
/** pid of particle 2 */
extern int    meta_pid2;
/** bias height */
extern double meta_bias_height;
/** bias width */
extern double meta_bias_width;

/** REACTION COORDINATE */
/** RC min */
extern double meta_xi_min;
/** RC max */
extern double meta_xi_max;
/** Force at boundaries */
extern double meta_f_bound;
/** Number of bins of RC */
extern int    meta_xi_num_bins;
/** Step between two bins of RC (determined by meta_xi_num_bins */
extern double meta_xi_step;


/* Arrays will be defined in metadynamics_init() */
/** Accumulated force array */
extern double *meta_acc_force;
/** Accumulated free energy profile */
extern double *meta_acc_fprofile;

/* Current vector of the reaction coordinate */
extern double *meta_cur_xi;
/* Current value of the reaction coordinate (scalar) */
extern double meta_val_xi;
/* Direction of the force that is applied (normalized) */
extern double *meta_apply_direction;

/*********************************
* functions
*********************************/

/** Initialize metadynamics on start of integration 
*  Create arrays if necessary. */
void meta_init();
/** Metadynamics main function:
* - Calculate reaction coordinate 
* - Update profile and biased force 
* - apply external force
*/
void meta_perform();

/** Calculate Lucy function */
double calculate_lucy(double xi, double xi_0);
/** Calculate derivative of Lucy function */
double calculate_deriv_lucy(double xi, double xi_0);

#endif

#endif
