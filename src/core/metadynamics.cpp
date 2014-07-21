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

#include "metadynamics.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "cells.hpp"

/** \file metadynamics.hpp 
*
*  This file contains routines to perform metadynamics.  Right now, the
*  reaction coordinate is defined between two particles (either distance 
*  or z-projected distance). Note that these
*  particles can be virtual sites, in order to handle molecules.
*
*  - set metadynamics options 
*  - initialize bias forces and free energy profiles 
*  - calculate reaction coordinate for each integration step 
*  - apply bias force on particles
*/


#ifdef METADYNAMICS
/* metadynamics switch */
int    meta_switch       = META_OFF;
/** pid of particle 1 */
int    meta_pid1         =       -1;
/** pid of particle 2 */
int    meta_pid2         =       -1;
/** bias height */
double meta_bias_height  =    0.001;
/** bias width */
double meta_bias_width   =      0.5;

/** REACTION COORDINATE */
/** RC min */
double meta_xi_min       =        1;
/** RC max */
double meta_xi_max       =        0;
/** Force at boundaries */
double meta_f_bound      =       10;
/** Number of bins of RC */
int    meta_xi_num_bins  =      100;
double meta_xi_step      =        1;


/** Accumulated force array */
double *meta_acc_force   =     NULL;
/** Accumulated free energy profile */
double *meta_acc_fprofile=     NULL;


double *meta_cur_xi      =     NULL;
double meta_val_xi       =       0.;
double *meta_apply_direction = NULL;

void meta_init(){
   if(meta_switch == META_OFF)   return;
   
   
   /* Initialize arrays if they're empty. These get freed upon calling the Tcl
    * parser */
   if (meta_acc_force == NULL || meta_acc_fprofile == NULL) {
      meta_acc_force       = (double *) calloc(meta_xi_num_bins * sizeof *meta_acc_force, sizeof *meta_acc_force);
      meta_acc_fprofile    = (double *) calloc(meta_xi_num_bins * sizeof *meta_acc_fprofile, sizeof *meta_acc_fprofile);
      meta_cur_xi          = (double *) calloc(3 * sizeof *meta_cur_xi, sizeof *meta_cur_xi);
      meta_apply_direction = (double *) calloc(3 * sizeof *meta_apply_direction, sizeof *meta_apply_direction);   
   }
   
   /* Check that the simulation uses onle a single processor. Otherwise exit. 
   *  MPI interface *not* implemented. */
   if (n_nodes != 1) {
      char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"Can't use metadynamics on more than one processor.\n");
      return;
   }
   
   meta_xi_step = (meta_xi_max-meta_xi_min)/(1.*meta_xi_num_bins);
}


/** Metadynamics main function:
* - Calculate reaction coordinate
* - Update profile and biased force
* - apply external force
*/
void meta_perform()
{
   if(meta_switch == META_OFF)  return;
   
   double ppos1[3] = {0,0,0}, ppos2[3] = {0,0,0}, factor;
   int img1[3], img2[3], c, i, np, flag1 = 0, flag2 = 0;
   Particle *p, *p1 = NULL, *p2 = NULL;
   Cell *cell;
   
   for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
         if (p[i].p.identity == meta_pid1) {
            flag1 = 1;
            p1 = &p[i];
            memcpy(ppos1, p[i].r.p, 3*sizeof(double));
            memcpy(img1, p[i].l.i, 3*sizeof(int));
            unfold_position(ppos1, img1);
            if (flag1 && flag2) {
               /* vector r2-r1 - Not a minimal image! Unfolded position */
               vector_subt(meta_cur_xi,ppos2,ppos1);
               break;
            }
         }
         if (p[i].p.identity == meta_pid2) {
            flag2 = 1;
            p2 = &p[i];
            memcpy(ppos2, p[i].r.p, 3*sizeof(double));
            memcpy(img2, p[i].l.i, 3*sizeof(int));
            unfold_position(ppos2, img2);
            if (flag1 && flag2) {
               /* vector r2-r1 - Not a minimal image! Unfolded position */
               vector_subt(meta_cur_xi,ppos2,ppos1);
               break;
            }
         }
      }
   }
   
   if (flag1 == 0 || flag2 == 0) {
      char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"Metadynamics: can't find pid1 or pid2.\n");
      return;
   }
   
   /* Now update free energy profile 
   * Here, we're following the functional form of 
   * Marsili etal., J Comp. Chem, 31 (2009).
   * Instead of gaussians, we use so-called Lucy's functions */
   
   for (i = 0; i < meta_xi_num_bins; ++i) {
      if (meta_switch == META_DIST) {
         // reaction coordinate value
         meta_val_xi = sqrt(sqrlen(meta_cur_xi));
         // Update free energy profile and biased force
         meta_acc_fprofile[i] -= calculate_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);
         meta_acc_force[i] -= calculate_deriv_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);
         
         // direction of the bias force
         unit_vector(meta_cur_xi,meta_apply_direction);
      } else if (meta_switch == META_REL_Z) {
         // reaction coordinate value: relative height of z_pid1 with respect to z_pid2
         meta_val_xi = -1.*meta_cur_xi[2];
         // Update free energy profile and biased force
         meta_acc_fprofile[i] -= calculate_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);
         meta_acc_force[i] -= calculate_deriv_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);
         
         // direction of the bias force (-1 to be consistent with META_DIST: from 1 to 2)
         meta_apply_direction[0] = meta_apply_direction[1] = 0.;
         meta_apply_direction[2] = -1.;
      } else {
         char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
         ERROR_SPRINTF(errtxt,"Undefined metadynamics scheme.\n");
         return;      
      }
   }
   
   /** Apply force */
   
   // Calculate the strength of the applied force
   if (meta_val_xi < meta_xi_min) {
      // below the lower bound
      factor = -1. * meta_f_bound * (meta_xi_min-meta_val_xi)/meta_xi_step;
   } else if (meta_val_xi > meta_xi_max) {
      // above the upper bound
      factor = meta_f_bound * (meta_val_xi-meta_xi_max)/meta_xi_step;
   } else {
      // within the RC interval
      i = (int) dround((meta_val_xi-meta_xi_min)/meta_xi_step);
      if (i < 0) i = 0;
      if (i >= meta_xi_num_bins) i=meta_xi_num_bins-1;
      factor = meta_acc_force[i];
   }
   
   /* cancel previous force to external force of particle */
   for (i = 0; i<3; ++i) {
      p1->f.f[i] +=       factor * meta_apply_direction[i];
      p2->f.f[i] += -1. * factor * meta_apply_direction[i];
   }
}


/** Calculate Lucy's function */
double calculate_lucy(double xi, double xi_0)
{  
   double dist = fabs(xi-xi_0);
   if (dist <= meta_bias_width) {
      return meta_bias_height 
      * (1+2*dist/meta_bias_width) 
      * pow(1-dist/meta_bias_width,2);
   } else 
      return 0.;
}

/** Calculate derivative of Lucy function */
double calculate_deriv_lucy(double xi, double xi_0)
{
   double dist = fabs(xi-xi_0);
   if (dist <= meta_bias_width) {
      double result = -2*meta_bias_height/meta_bias_width
      * (pow(1-dist/meta_bias_width,2)
      -(1+2*dist/meta_bias_width)*(1-dist/meta_bias_width));
      if (xi < xi_0)
         result *= -1.;
      return result;
   } else
      return 0.;
}



#endif
