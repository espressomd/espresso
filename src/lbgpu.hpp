/* 
   Copyright (C) 2010,2011,2012,2013 The ESPResSo project

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

/** \file lbgpu.hpp
 * Header file for lbgpu.cpp
 *
 * This is the header file for the Lattice Boltzmann implementation in lbgpu_cfile.cpp
 */

#ifndef LB_GPU_H
#define LB_GPU_H

#include "utils.hpp"
#include "config.hpp"

#ifdef LB_GPU

/* For the D3Q19 model most functions have a separate implementation
 * where the coefficients and the velocity vectors are hardcoded
 * explicitly. This saves a lot of multiplications with 1's and 0's
 * thus making the code more efficient. */
#define D3Q19
#define LBQ 19

/** Note these are usef for binary logic so should be powers of 2 */
#define LB_COUPLE_NULL        1
#define LB_COUPLE_TWO_POINT   2
#define LB_COUPLE_THREE_POINT 4

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
#ifdef SHANCHEN
#define LBPAR_COUPLING 8
#define LBPAR_MOBILITY 9
#endif
/*@}*/

/**-------------------------------------------------------------------------*/
/** Data structure holding the parameters for the Lattice Boltzmann system for gpu. */
typedef struct {
  /** number density (LJ units) */
  float rho[LB_COMPONENTS];
  /** mu (LJ units) */
  float mu[LB_COMPONENTS];
  /*viscosity (LJ) units */
  float viscosity[LB_COMPONENTS];
  /** relaxation rate of shear modes */
  float gamma_shear[LB_COMPONENTS];
  /** relaxation rate of bulk modes */
  float gamma_bulk[LB_COMPONENTS];
  /**      */
  float gamma_odd[LB_COMPONENTS];
  float gamma_even[LB_COMPONENTS];
  /** friction coefficient for viscous coupling (LJ units)
   * Note that the friction coefficient is quite high and may
   * lead to numerical artifacts with low order integrators */
  float friction[LB_COMPONENTS];
  /** amplitude of the fluctuations in the viscous coupling */
  /** Switch indicating what type of coupling is used, can either
  use nearest neighbors or next nearest neighbors. */
  int lb_couple_switch;

  float lb_coupl_pref[LB_COMPONENTS];
  float lb_coupl_pref2[LB_COMPONENTS];
  float bulk_viscosity[LB_COMPONENTS];

  /** lattice spacing (LJ units) */
  float agrid;

  /** time step for fluid propagation (LJ units)
   *  Note: Has to be larger than MD time step! */
  float tau;

  /** MD timestep */
  float time_step;

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

  unsigned int reinit;

#ifdef SHANCHEN
  /** mobility. They are actually LB_COMPONENTS-1 in number, we leave LB_COMPONENTS here for practical reasons*/
  float gamma_mobility[LB_COMPONENTS];
  float mobility[LB_COMPONENTS];
  float coupling[LB_COMPONENTS*LB_COMPONENTS];
  int remove_momentum;
#endif // SHANCHEN  

} LB_parameters_gpu;

/** Data structure holding the conserved quantities for the Lattice Boltzmann system. */
typedef struct {

  /** density of the node */
  float rho[LB_COMPONENTS];
  /** veolcity of the node */

  float v[3];

} LB_rho_v_gpu;
/* this structure is almost duplicated for memory efficiency. When the stress 
   tensor element are needed at every timestep, this features should be explicitly
   switched on */
typedef struct { 
  /** density of the node */
  float rho[LB_COMPONENTS];
  /** veolcity of the node */
  float v[3];
  /** pressure tensor */
  float pi[6];  
} LB_rho_v_pi_gpu;

/** Data structure holding the velocity densities for the Lattice Boltzmann system. */
typedef struct {

  /** velocity density of the node */
  float *vd;
  /** seed for the random gen */
  unsigned int *seed;
  /** flag indicating whether this site belongs to a boundary */
  unsigned int *boundary;

} LB_nodes_gpu;

/** Data structure for the randomnr and the seed. */
typedef struct {

  float randomnr[2];

  unsigned int seed;

} LB_randomnr_gpu;

typedef struct {

  float *force;

} LB_node_force_gpu;

typedef struct {

  float force[3];

  unsigned int index;

} LB_extern_nodeforce_gpu;


void on_lb_params_change_gpu(int field);

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** 
 */

/** Switch indicating momentum exchange between particles and fluid */
extern LB_parameters_gpu lbpar_gpu;
extern LB_rho_v_pi_gpu *host_values;
extern int transfer_momentum_gpu;
extern LB_extern_nodeforce_gpu *extern_nodeforces_gpu;
extern int n_lb_boundaries;
#ifdef ELECTROKINETICS
extern LB_node_force_gpu node_f;
extern int ek_initialized;
#endif


/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/


void lb_get_device_values_pointer(LB_rho_v_gpu** pointeradress);
void lb_get_boundary_force_pointer(float** pointeradress);
void lb_get_lbpar_pointer(LB_parameters_gpu** pointeradress);
void lb_get_para_pointer(LB_parameters_gpu** pointeradress);
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
void lb_realloc_particles_GPU_leftovers(LB_parameters_gpu *lbpar_gpu);

void lb_init_GPU(LB_parameters_gpu *lbpar_gpu);
void lb_integrate_GPU();
#ifdef SHANCHEN
void lb_calc_shanchen_GPU();
void lattice_boltzmann_calc_shanchen_gpu();
#endif
void lb_free_GPU();
void lb_get_values_GPU(LB_rho_v_pi_gpu *host_values);
void lb_print_node_GPU(int single_nodeindex, LB_rho_v_pi_gpu *host_print_values);
#ifdef LB_BOUNDARIES_GPU
void lb_init_boundaries_GPU(int n_lb_boundaries, int number_of_boundnodes, int* host_boundary_node_list, int* host_boundary_index_list, float* lb_bounday_velocity);
#endif
void lb_init_extern_nodeforces_GPU(int n_extern_nodeforces, LB_extern_nodeforce_gpu *host_extern_nodeforces, LB_parameters_gpu *lbpar_gpu);

void lb_calc_particle_lattice_ia_gpu();

void lb_calc_fluid_mass_GPU(double* mass);
void lb_calc_fluid_momentum_GPU(double* host_mom);
void lb_remove_fluid_momentum_GPU(void);
void lb_calc_fluid_temperature_GPU(double* host_temp);
void lb_get_boundary_flag_GPU(int single_nodeindex, unsigned int* host_flag);
void lb_get_boundary_flags_GPU(unsigned int* host_bound_array);

void lb_set_node_velocity_GPU(int single_nodeindex, float* host_velocity);
void lb_set_node_rho_GPU(int single_nodeindex, float* host_rho);

void reinit_parameters_GPU(LB_parameters_gpu *lbpar_gpu);
void lb_reinit_extern_nodeforce_GPU(LB_parameters_gpu *lbpar_gpu);
void lb_reinit_GPU(LB_parameters_gpu *lbpar_gpu);
int lb_lbnode_set_extforce_GPU(int ind[3], double f[3]);
void lb_gpu_get_boundary_forces(double* forces);
void lb_save_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, unsigned int *host_checkpoint_boundary, float *host_checkpoint_force);
void lb_load_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, unsigned int *host_checkpoint_boundary, float *host_checkpoint_force);


/*@{*/

#endif /* LB_GPU */

#endif /* LB_GPU_H */
