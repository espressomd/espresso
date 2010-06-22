// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#ifndef METADYNAMICS_H
#define METADYNAMICS_H

#include <string.h>
#include <tcl.h>
#include <math.h>
#include "particle_data.h"
#include "utils.h"
#include "parser.h"
#include "communication.h"
#include "cells.h"

/** \file metadynamics.h 
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

/** Implementation of the Tcl command \ref tcl_metadynamics. This function
 *  allows to change the parameters of metadynamics */
int metadynamics(ClientData data, Tcl_Interp *interp, int argc, char **argv);
/** Print metadynamics options and parameters */
int meta_print(Tcl_Interp *interp);
int meta_usage(Tcl_Interp *interp, int argc, char **argv);
int meta_parse_off(Tcl_Interp *interp, int argc, char **argv);
/** Reaction coordinates available */
int meta_parse_distance(Tcl_Interp *interp, int argc, char **argv);
int meta_parse_relative_z(Tcl_Interp *interp, int argc, char **argv);
/** Input/Output stuff */
int meta_print_stat(Tcl_Interp *interp, int argc, char **argv);
int meta_parse_stat(Tcl_Interp *interp, int argc, char **argv);

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
