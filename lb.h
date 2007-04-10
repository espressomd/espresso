/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
 */

/** \file lb.h
 * Header file for lb.c
 *
 * This is the header file for the Lattice Boltzmann implementation in lb.c
 */

#ifndef LB_H
#define LB_H

#include <tcl.h>
#include "utils.h"
#include "lattice.h"

#ifdef LB

/* For the D3Q19 model most functions have a separate implementation
 * where the coefficients and the velocity vectors are hardcoded
 * explicitly. This saves a lot of multiplications with 1's and 0's
 * thus making the code more efficient. */
#define D3Q19

/** \name Parameter fields for Lattice Boltzmann 
 * The numbers are referenced in \ref mpi_bcast_lb_params 
 * to determine what actions have to take place upon change
 * of the respective parameter. */
/*@{*/
#define LBPAR_DENSITY   0 /**< fluid density */
#define LBPAR_VISCOSITY 1 /**< fluid kinematic viscosity */
#define LBPAR_AGRID     2 /**< grid constant for fluid lattice */
#define LBPAR_TAU       3 /**< time step for fluid propagation */
#define LBPAR_FRICTION  4 /**< friction coefficient for viscous coupling between particles and fluid */
#define LBPAR_EXTFORCE  5 /**< external force acting on the fluid */
/*@}*/

/** Description of the LB Model in terms of the unit vectors of the 
 *  velocity sub-lattice and the corresponding coefficients 
 *  of the pseudo-equilibrium distribution */
typedef struct {

  /** number of velocities */
  int n_veloc ;

  /** unit vectors of the velocity sublattice */
  double (*c)[3] ;

  /** coefficients in the pseudo-equilibrium distribution */
  double (*coeff)[4];

  /** weights in the functional for the equilibrium distribution */
  double (*w);

  /** speed of sound squared */
  double c_sound_sq;

} LB_Model;

/** Data structure for fluid on a local lattice site */
typedef struct {

  /** local populations of the velocity directions */
  double *n;

  /** local density */
  double *rho;

  /** local momentum */
  double *j;

  /** local stress tensor */
  double *pi;

#ifdef CONSTRAINTS
   /** flag indicating whether this site belongs to a boundary */
   int boundary;

  /** normal vector of the boundary surface */
  double *nvec;
#endif

#ifndef D3Q19
  /** temporary storage for the new populations */
  double *n_new;
#endif

} LB_FluidNode;

/** Data structure holding the parameters for the Lattice Boltzmann system. */
typedef struct {

  /** number density (LJ units) */
  double rho;

  /** kinematic viscosity (LJ units) */
  double viscosity;

  /** bulk viscosity (LJ units) */
  double bulk_viscosity;

  /** lattice spacing (LJ units) */
  double agrid;

  /** time step for fluid propagation (LJ units)
   *  Note: Has to be larger than MD time step! */
  double tau;

  /** friction coefficient for viscous coupling (LJ units)
   * Note that the friction coefficient is quite high and may
   * lead to numerical artifacts with low order integrators */
  double friction;

  /** external force applied to the fluid at each lattice site (LJ units) */
  double ext_force[3];
          
} LB_Parameters;

/** The DnQm model to be used. */
extern LB_Model lbmodel;

/** Struct holding the Lattice Boltzmann parameters */
extern LB_Parameters lbpar; 

/** The underlying lattice */
extern Lattice lblattice;

/** Pointer to the fluid nodes */
extern LB_FluidNode *lbfluid;

/** Switch indicating momentum exchange between particles and fluid */
extern int transfer_momentum;

extern double lblambda;

extern double lblambda_bulk;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init();

/** (Re-)initializes the derived parameters
 *  for the Lattice Boltzmann system.
 *  The current state of the fluid is unchanged. */
void lb_reinit_parameters();

/** (Re-)initializes the fluid. */
void lb_reinit_fluid();

/** Sets the density and momentum on a local lattice site.
 * @param index The index of the lattice site within the local domain (Input)
 * @param rho   Local density of the fluid (Input)
 * @param j     Local momentum of the fluid (Input)
 */
void lb_set_local_fields(int index, const double rho, const double *v);

/** Returns the mass, momentum and stress of a local lattice site.
 * @param index The index of the lattice site within the local domain (Input)
 * @param rho   Local density of the fluid (Output)
 * @param j     Local momentum of the fluid (Output)
 * @param pi    Local stress tensor of the fluid (Output)
 */
void lb_get_local_fields(int index, double *rho, double *j, double *pi);

/** Propagates the Lattice Boltzmann system for one time step.
 * This function performs the collision step and the streaming step.
 * If external forces are present, they are applied prior to the collisions.
 * If boundaries are present, it also applies the boundary conditions.
 */
void lb_propagate();

/** Calculates the coupling of MD particles to the LB fluid.
 * This function  is called from \ref force_calc. The force is added
 * to the particle force and the corresponding momentum exchange is
 * applied to the fluid. 
 * Note that this function changes the state of the fluid!
 */
void calc_particle_lattice_ia();

/** Calculate the local fluid density.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param local_node The local lattice site (Input).
 */
MDINLINE void lb_calc_local_rho(LB_FluidNode *local_node) {

  double *local_n = local_node->n;
  double *local_rho = local_node->rho;
  double avg_rho = lbpar.rho*lbpar.agrid*lbpar.agrid*lbpar.agrid;

#ifdef D3Q19
  *local_rho =   avg_rho
               + local_n[0] 
               + local_n[1]  + local_n[2]  
               + local_n[3]  + local_n[4]  
               + local_n[5]  + local_n[6] 
               + local_n[7]  + local_n[8]  + local_n[9]  + local_n[10] 
               + local_n[11] + local_n[12] + local_n[13] + local_n[14] 
               + local_n[15] + local_n[16] + local_n[17] + local_n[18];
#else
  int i;
  *local_rho = 0.0;
  for (i=0;i<n_veloc;i++) {
    *local_rho += local_n[i] + c[i][0]*avg_rho;
  }
#endif

}

/** Calculate the local fluid momentum.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param local_node The local lattice site (Input).
 */
MDINLINE void lb_calc_local_j(LB_FluidNode *local_node) {

  double *local_n = local_node->n;
  double *local_j = local_node->j;

#ifdef D3Q19
  local_j[0] =   local_n[1]  - local_n[2] 
               + local_n[7]  - local_n[8]  + local_n[9]  - local_n[10] 
               + local_n[11] - local_n[12] + local_n[13] - local_n[14];
  local_j[1] =   local_n[3]  - local_n[4]
               + local_n[7]  - local_n[8]  - local_n[9]  + local_n[10]
               + local_n[15] - local_n[16] + local_n[17] - local_n[18]; 
  local_j[2] =   local_n[5]  - local_n[6]  
               + local_n[11] - local_n[12] - local_n[13] + local_n[14]
               + local_n[15] - local_n[16] - local_n[17] + local_n[18];
#else
  int i;
  double tmp;
  local_j[0] = 0.0;
  local_j[1] = 0.0;
  local_j[2] = 0.0;
  for (i=0;i<n_veloc;i++) {
    tmp = local_n[i] + c[i][0]*avg_rho;
    local_j[0] += c[i][0] * tmp;
    local_j[1] += c[i][1] * tmp;
    local_j[2] += c[i][2] * tmp;
  }
#endif

}

/** Calculate the local fluid stress.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param local_node The local lattice site (Input).
 */
MDINLINE void lb_calc_local_pi(LB_FluidNode *local_node) {

  double *local_n  = local_node->n;
  double *local_pi = local_node->pi;
  double avg_rho = lbpar.rho*lbpar.agrid*lbpar.agrid*lbpar.agrid;
    
#ifdef D3Q19
  local_pi[0] =   avg_rho/3.0
                + local_n[1]  + local_n[2]  
                + local_n[7]  + local_n[8]  + local_n[9]  + local_n[10] 
                + local_n[11] + local_n[12] + local_n[13] + local_n[14];
  local_pi[2] =   avg_rho/3.0
                + local_n[3]  + local_n[4]  
                + local_n[7]  + local_n[8]  + local_n[9]  + local_n[10]
                + local_n[15] + local_n[16] + local_n[17] + local_n[18];
  local_pi[5] =   avg_rho/3.0
                + local_n[5]  + local_n[6]  
                + local_n[11] + local_n[12] + local_n[13] + local_n[14] 
                + local_n[15] + local_n[16] + local_n[17] + local_n[18];
  local_pi[1] =   local_n[7]  + local_n[8]  - local_n[9]  - local_n[10];
  local_pi[3] =   local_n[11] + local_n[12] - local_n[13] - local_n[14];
  local_pi[4] =   local_n[15] + local_n[16] - local_n[17] - local_n[18];
#else
  int i;
  double tmp;
  local_pi[0] = 0.0;
  local_pi[1] = 0.0;
local_pi[2] = 0.0;
  local_pi[3] = 0.0;
  local_pi[4] = 0.0;
  local_pi[5] = 0.0;
  for (i=0;i<n_veloc;i++) {
    tmp = local_n[i] + c[i][0]*avg_rho;
    local_pi[0] += c[i][0] * c[i][0] * tmp;
    local_pi[1] += c[i][0] * c[i][1] * tmp;
    local_pi[2] += c[i][1] * c[i][1] * tmp;
    local_pi[3] += c[i][0] * c[i][2] * tmp;
    local_pi[4] += c[i][1] * c[i][2] * tmp;
    local_pi[5] += c[i][2] * c[i][2] * tmp;
  }
#endif

}

/** Calculate the local fluid fields.
 * The calculation is implemented explicitly for the special case of D3Q19.
 *
 * Original Author: Ahlrichs 06/11/97, 29/03/98
 *
 * @param local_node   The local lattice site.
 * @param calc_pi_flag Flag indicating whether stress tensor should be
 *                     computed.
 */
MDINLINE void lb_calc_local_fields(LB_FluidNode *local_node,int calc_pi_flag) {

  double *local_n   = local_node->n;
  double *local_rho = local_node->rho;
  double *local_j   = local_node->j;
  double *local_pi  = local_node->pi;
  double avg_rho = lbpar.rho*lbpar.agrid*lbpar.agrid*lbpar.agrid;

#ifdef D3Q19
  *local_rho =   avg_rho
               + local_n[0]  
               + local_n[1]  + local_n[2]  
               + local_n[3]  + local_n[4] 
               + local_n[5]  + local_n[6]  
               + local_n[7]  + local_n[8]  + local_n[9]  + local_n[10] 
               + local_n[11] + local_n[12] + local_n[13] + local_n[14]
               + local_n[15] + local_n[16] + local_n[17] + local_n[18];

  local_j[0] =   local_n[1]  - local_n[2]
               + local_n[7]  - local_n[8]  + local_n[9]  - local_n[10]
               + local_n[11] - local_n[12] + local_n[13] - local_n[14];
  local_j[1] =   local_n[3]  - local_n[4]
               + local_n[7]  - local_n[8]  - local_n[9]  + local_n[10]
               + local_n[15] - local_n[16] + local_n[17] - local_n[18]; 
  local_j[2] =   local_n[5]  - local_n[6]
               + local_n[11] - local_n[12] - local_n[13] + local_n[14]
               + local_n[15] - local_n[16] - local_n[17] + local_n[18];
  
  if (calc_pi_flag) {
    local_pi[0] =   avg_rho/3.0
                  + local_n[1]  + local_n[2]  
                  + local_n[7]  + local_n[8]  + local_n[9]  + local_n[10]
                  + local_n[11] + local_n[12] + local_n[13] + local_n[14];
    local_pi[2] =   avg_rho/3.0
                  + local_n[3]  + local_n[4]
                  + local_n[7]  + local_n[8]  + local_n[9]  + local_n[10]
                  + local_n[15] + local_n[16] + local_n[17] + local_n[18];
    local_pi[5] =   avg_rho/3.0
                  + local_n[5]  + local_n[6]
                  + local_n[11] + local_n[12] + local_n[13] + local_n[14]
                  + local_n[15] + local_n[16] + local_n[17] + local_n[18];
    local_pi[1] =   local_n[7]  - local_n[9] + local_n[8] - local_n[10];
    local_pi[3] =   local_n[11] + local_n[12] - local_n[13] - local_n[14];
    local_pi[4] =   local_n[15] + local_n[16] - local_n[17] - local_n[18];

  }
#else
  int i;
  double nlocal;

  *local_rho = 0.0;

  local_j[0] = 0.0;
  local_j[1] = 0.0;
  local_j[2] = 0.0;

  if (calc_pi_flag) {
    local_pi[0] = 0.0;
    local_pi[1] = 0.0;
    local_pi[2] = 0.0;
    local_pi[3] = 0.0;
    local_pi[4] = 0.0;
    local_pi[5] = 0.0;
  }

  for (i=0;i<n_veloc;i++) {
    tmp = local_n[i] + c[i][0];
    
    *local_rho += tmp;

    local_j[0] += c[i][0] * tmp;
    local_j[1] += c[i][1] * tmp;
    local_j[2] += c[i][2] * tmp;

    if (calc_pi_flag) {
      local_pi[0] += c[i][0] * c[i][0] * tmp;
      local_pi[1] += c[i][0] * c[i][1] * tmp;
      local_pi[2] += c[i][1] * c[i][1] * tmp;
      local_pi[3] += c[i][0] * c[i][2] * tmp;
      local_pi[4] += c[i][1] * c[i][2] * tmp;
      local_pi[5] += c[i][2] * c[i][2] * tmp;
    }

  }
#endif

}

#endif /* LB */

/** Parser for the TCL command \ref lbfluid. */
int lbfluid_cmd(ClientData data, Tcl_Interp *interp, int argc, char **argv) ;

#endif /* LB_H */

/*@}*/
