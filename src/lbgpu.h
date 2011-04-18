/* $Id: lbgpu.h $
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

/** \file lbgpu.h
 * Header file for lbgpu.c
 *
 * This is the header file for the Lattice Boltzmann implementation in lbgpu_cfile.c
 */




#ifndef LB_GPU_H
#define LB_GPU_H

#include <tcl.h>
#include "utils.h"
//#include "lattice.h"

#ifdef LB_GPU

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
#define LBPAR_BULKVISC  6 /**< fluid bulk viscosity */
#ifdef CONSTRAINTS
#define LBPAR_BOUNDARY  7 /**< boundary parameters */
#endif
/*@}*/


/**-------------------------------------------------------------------------*/
/** Data structure holding the parameters for the Lattice Boltzmann system for gpu. */
typedef struct {

  /** number density (LJ units) */
  float rho;

  /** mu (LJ units) */
  float mu;

  /*viscosity (LJ) units */
  float viscosity;

  /** relaxation rate of shear modes */
  float gamma_shear;
  /** relaxation rate of bulk modes */
  float gamma_bulk;
  /**      */
  float gamma_odd;
  float gamma_even;

  /** lattice spacing (LJ units) */
  float agrid;

  /** time step for fluid propagation (LJ units)
   *  Note: Has to be larger than MD time step! */
  float tau;

  /** friction coefficient for viscous coupling (LJ units)
   * Note that the friction coefficient is quite high and may
   * lead to numerical artifacts with low order integrators */
  float friction;
  /** MD tiemstep */
  float time_step;
  /** amplitude of the fluctuations in the viscous coupling */
  float lb_coupl_pref;

  float lb_coupl_pref2;

  float bulk_viscosity;

  unsigned int dim_x;
  unsigned int dim_y;
  unsigned int dim_z;

  unsigned int number_of_nodes;
  unsigned int number_of_particles;
  /** Flag indicating whether fluctuations are present. */
  int fluct;
  /**to calc and print out phys values */
  int calc_val;

  int external_force;

  float ext_force[3];

  unsigned int your_seed;

} LB_parameters_gpu;

/** Data structure holding the velocitydensities for the Lattice Boltzmann system. */
typedef struct {

  /** velocitydensity of the node */
  float *vd;
  /** seed for the random gen */
  unsigned int *seed;
  /** flag indicating whether this site belongs to a boundary */
  unsigned int *boundary;

} LB_nodes_gpu;

/** Data structure holding the phys. values for the Lattice Boltzmann system. */
typedef struct {

  /** velocitydensity of the node */
  float rho;

  /** veolcity of the node */
  float v[3];

  /** stresstensor of the node */
  /** use this value only (due to memory saving) if you want to print out the value (used in calc_values)*/
  //float pi[6];

} LB_values_gpu;

/** Data structure for the randomnr and the seed. */
typedef struct {

  float randomnr[2];

  unsigned int seed;

} LB_randomnr_gpu;

typedef struct {
  /** force on the particle given to md part */
  float f[3];

} LB_particle_force_gpu;

typedef struct {
  /** particle position given from md part*/
  float p[3];
  /** particle momentum struct velocity p.m->v*/
  float v[3];
#ifdef LB_ELECTROHYDRODYNAMICS
  float mu_E[2];
#endif
  unsigned int fixed;

} LB_particle_gpu;

typedef struct {

  float *force;

} LB_node_force_gpu;

typedef struct {

  float force[3];

  unsigned int index;

} LB_extern_nodeforce_gpu;

typedef struct {

  unsigned int seed;

} LB_particle_seed_gpu;

#ifdef __cplusplus
extern "C" {
#endif

extern LB_parameters_gpu lbpar_gpu;

/** Switch indicating momentum exchange between particles and fluid */
extern int transfer_momentum_gpu;

extern unsigned int lb_boundaries_bb_gpu;

extern LB_extern_nodeforce_gpu *extern_nodeforces_gpu;

#ifdef __cplusplus
}
#endif
/** Eigenvalue of collision operator corresponding to shear viscosity. */
//extern double lblambda;

/** Eigenvalue of collision operator corresponding to bulk viscosity. */
//extern double lblambda_bulk;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Updates the Lattice Boltzmann system for one time step.
 * This function performs the collision step and the streaming step.
 * If external forces are present, they are applied prior to the collisions.
 * If boundaries are present, it also applies the boundary conditions.
 */
#ifdef __cplusplus
extern "C" {
#endif
void lattice_boltzmann_update_gpu();

/** (Pre-)initializes data structures. */
void lb_pre_init_gpu();

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init_gpu();

/** (Re-)initializes the derived parameters
 *  for the Lattice Boltzmann system.
 *  The current state of the fluid is unchanged. */
void lb_reinit_parameters_gpu();

/** (Re-)initializes the fluid. */
void lb_reinit_fluid_gpu();

/** Resets the forces on the fluid nodes */
//void lb_reinit_forces();

/** (Re-)initializes the particle array*/
void lb_realloc_particles_gpu();


void lb_init_GPU(LB_parameters_gpu *lbpar_gpu);
void lb_integrate_GPU();
void lb_particle_GPU(LB_particle_gpu *host_data);
void lb_free_GPU();
void lb_get_values_GPU(LB_values_gpu *host_values);
void lb_realloc_particle_GPU(LB_parameters_gpu *lbpar_gpu, LB_particle_gpu **host_data);
void lb_copy_forces_GPU(LB_particle_force_gpu *host_forces);
void lb_print_node_GPU(int single_nodeindex, LB_values_gpu *host_print_values);
void lb_init_boundaries_GPU(int number_of_boundnodes, int *host_boundindex);
void lb_init_extern_nodeforces_GPU(int n_extern_nodeforces, LB_extern_nodeforce_gpu *host_extern_nodeforces, LB_parameters_gpu *lbpar_gpu);

void lb_calc_particle_lattice_ia_gpu();
void lb_send_forces_gpu();
void calc_fluid_momentum_GPU(double* mom);
void calc_fluid_temperature_GPU(double* cpu_temp);

#ifdef __cplusplus
}
#endif

/** Parser for the TCL command lbnode. */
int tclcommand_lbnode_gpu(Tcl_Interp *interp, int argc, char **argv);

/** Parser for the TCL command \ref lbfluid. */
int tclcommand_lbfluid_gpu(Tcl_Interp *interp, int argc, char **argv);

int tclcommand_lbnode_extforce_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv);

int tclcommand_lbprint_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv);
#endif /* LB_GPU */
#endif /* LB_GPU_H */

/*@}*/
