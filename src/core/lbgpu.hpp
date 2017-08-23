/* 
   Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

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
#include "observables/profiles.hpp" 
#include "particle_data.hpp"
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

#if defined(LB_DOUBLE_PREC) || defined(EK_DOUBLE_PREC)
typedef double lbForceFloat;
#else
typedef float lbForceFloat;
#endif
/**-------------------------------------------------------------------------*/
/** Data structure holding the parameters for the Lattice Boltzmann system for gpu. */
typedef struct {
  /** number density (LJ units) */
  float rho[LB_COMPONENTS];
  /** mu (LJ units) */
  float mu[LB_COMPONENTS];
  /** viscosity (LJ) units */
  float viscosity[LB_COMPONENTS];
  /** relaxation rate of shear modes */
  float gamma_shear[LB_COMPONENTS];
  /** relaxation rate of bulk modes */
  float gamma_bulk[LB_COMPONENTS];
  /**      */
  float gamma_odd[LB_COMPONENTS];
  float gamma_even[LB_COMPONENTS];
  /** flag determining whether gamma_shear, gamma_odd, and gamma_even are calculated
   *  from gamma_shear in such a way to yield a TRT LB with minimized slip at
   *  bounce-back boundaries */
  bool is_TRT;
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
#ifdef LB_BOUNDARIES_GPU
  unsigned int number_of_boundnodes;
#endif
  unsigned int number_of_anchors;
  /** Flag indicating whether fluctuations are present. */
  int fluct;
  /**to calc and print out phys values */
  int calc_val;

  int external_force;

  float ext_force[3*LB_COMPONENTS];

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
  /** Signed flag indicating whether this site belongs to a boundary. */
  int *boundary;
  /** Soon-to-be boundary flag needed for boundary propagation */
  int *boundary_buffer;
  
} LB_nodes_gpu;

 /**Data structure holding information necessary to integrate moving boundaries on gpu */
#endif // LB_GPU
typedef struct {
  float velocity[3];//stored as a/tau
  float omega[3];//stored as 1/tau, on host it's 1/1 
  float center[3];//folded coordinate. Absolute coord stored as double on host.

  //float orientation[3];
  float mass; //TOTAL mass (sum over anchor mass and central mass)
  float scaled_mass;//this is used for force calculations m/a^3
  float force_add[3];//momentum conservation, not used by default
  float force_hyd[3];//hydrodynamic interaction force
  float torque[3];//torque generated by hyd ia f
  float radius; //for sphere
  int n_anchors;
  std::vector<float> anchors; //stores sphere anchors rel. pos. and their radii (x,y,z,r) for arbitrary bounds
  int index; //which boundary this is (0 is fluid, then -1,-2,-3...)
  
 //The folowing values are only used to pass information to the host integrator on init
 //They need to be passed from parser via LB_Boundary struct to LB_moving_boundary struct
 //with set_moving_boundary_struct() in lb-boundaries.cpp
 //then set to a host Particle struct with set_virtual_...() in lbgpu_cuda.cu
 
  double quat[4]; //orientation is integrated on host only
#ifdef EXTERNAL_FORCES
  double ext_force[3];
  double body_force[3];
  double ext_torque[3];
  double body_torque[3];
#endif
  double rinertia[3];
} LB_moving_boundary;
#ifdef LB_GPU

/** Data structure for the randomnr and the seed. */
typedef struct {

  float randomnr[2];

  unsigned int seed;

} LB_randomnr_gpu;

typedef struct {

  lbForceFloat *force;
  float *scforce;
#if defined(IMMERSED_BOUNDARY) || defined(EK_DEBUG)

  // We need the node forces for the velocity interpolation at the virtual particles' position
  // However, LBM wants to reset them immediately after the LBM update
  // This variable keeps a backup
  lbForceFloat *force_buf;
#endif

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
#ifdef ELECTROKINETICS
extern LB_node_force_gpu node_f;
extern int ek_initialized;
#endif


/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

void lb_GPU_sanity_checks();

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
void reset_LB_forces_GPU(bool buffer = true);

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

/** Print moving boundaries */
void lb_convert_force_to_body_for_print(Particle *p, double *print_val);
void lb_convert_torque_to_space_for_print(Particle *p, double *print_val);
void lb_convert_omega_to_space_for_print(Particle *p, double *print_val);
void lb_convert_ext_torque_to_body_for_print(Particle *p, double *print_val);
int lb_lbfluid_print_moving_pos(int part_num, double *print_val);
int lb_lbfluid_print_moving_vel(int part_num, double *print_val);
int lb_lbfluid_print_moving_orientation(int part_num, double *print_val);
int lb_lbfluid_print_moving_omega_body(int part_num, double *print_val);
int lb_lbfluid_print_moving_omega_lab(int part_num, double *print_val);
int lb_lbfluid_print_moving_force_body(int part_num, double *print_val);
int lb_lbfluid_print_moving_force_lab(int part_num, double *print_val);
int lb_lbfluid_print_moving_torque_body(int part_num, double *print_val);
int lb_lbfluid_print_moving_torque_lab(int part_num, double *print_val);
int lb_lbfluid_print_moving_n(int print_val);


#ifdef LB_BOUNDARIES_GPU
void lb_init_boundaries_GPU(int n_lb_boundaries, int n_lb_moving_boundaries, int number_of_boundnodes, int* host_boundary_node_list, int* host_boundary_index_list, float* lb_bounday_velocity, LB_moving_boundary *host_moving_boundary);
void convert_omega_space_to_body(LB_moving_boundary *lbb, Particle *p); //rotation functions for init with either custom arguments or names
void convert_omega_anchors_body_to_space(Particle *p, float *omega, float *h_anchors);
void convert_force_body_to_space(Particle *p);
void lb_convert_quat_to_quatu(double quat[4], double quatu[3]);
#endif
void lb_init_extern_nodeforces_GPU(int n_extern_nodeforces, LB_extern_nodeforce_gpu *host_extern_nodeforces, LB_parameters_gpu *lbpar_gpu);

void lb_calc_particle_lattice_ia_gpu();

void lb_calc_fluid_mass_GPU(double* mass);
void lb_calc_fluid_momentum_GPU(double* host_mom);
void lb_remove_fluid_momentum_GPU(void);
void lb_calc_fluid_temperature_GPU(double* host_temp);
void lb_get_boundary_flag_GPU(int single_nodeindex, int* host_flag);
void lb_get_boundary_flags_GPU(int* host_bound_array);

void lb_set_node_velocity_GPU(int single_nodeindex, float* host_velocity);
void lb_set_node_rho_GPU(int single_nodeindex, float* host_rho);

void reinit_parameters_GPU(LB_parameters_gpu *lbpar_gpu);
void lb_reinit_extern_nodeforce_GPU(LB_parameters_gpu *lbpar_gpu);
void lb_reinit_GPU(LB_parameters_gpu *lbpar_gpu);
int lb_lbnode_set_extforce_GPU(int ind[3], double f[3]);
void lb_gpu_get_boundary_forces(double* forces);
void lb_save_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, int *host_checkpoint_boundary, lbForceFloat *host_checkpoint_force);
void lb_load_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, int *host_checkpoint_boundary, lbForceFloat *host_checkpoint_force);
int lb_lbfluid_save_checkpoint_wrapper(char* filename, int binary);
int lb_lbfluid_load_checkpoint_wrapper(char* filename, int binary);

void lb_lbfluid_remove_total_momentum();
void lb_lbfluid_fluid_add_momentum(float momentum[3]);
void lb_lbfluid_calc_linear_momentum(float momentum[3], int include_particles, int include_lbfluid);
void lb_lbfluid_particles_add_momentum(float velocity[3]);

void lb_lbfluid_set_population( int[3], float[LBQ], int );
void lb_lbfluid_get_population( int[3], float[LBQ], int );

//int statistics_observable_lbgpu_radial_velocity_profile(radial_profile_data* pdata, double* A, unsigned int n_A);
//int statistics_observable_lbgpu_velocity_profile(profile_data* pdata, double* A, unsigned int n_A);

/*@{*/

#endif /* LB_GPU */

#endif /* LB_GPU_H */
