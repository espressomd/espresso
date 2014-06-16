/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file lb.hpp
 * Header file for lb.cpp
 *
 * This is the header file for the Lattice Boltzmann implementation in lb.cpp
 */

#ifndef LB_H
#define LB_H

#include "utils.hpp"
#include "lattice.hpp"

extern int lb_components ; // global variable holding the number of fluid components

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
#define LBPAR_BULKVISC  6 /**< fluid bulk viscosity */

/** Note these are usef for binary logic so should be powers of 2 */
#define LB_COUPLE_NULL        1
#define LB_COUPLE_TWO_POINT   2
#define LB_COUPLE_THREE_POINT 4
  
/*@}*/
  /** Some general remarks:
   * This file implements the LB D3Q19 method to Espresso. The LB_Model
   * construction is preserved for historical reasons and might be removed
   * soon. It is constructed as a multi-relaxation time LB, thus all populations
   * are converted to modes, then collision is performed and transfered back
   * to population space, where the streaming is performed. 
   *
   * For performance reasons it is clever to do streaming and collision at the same time
   * because every fluid node has to be read and written only once. This increases
   * mainly cache efficiency. 
   * Two alternatives are implemented: stream_collide and collide_stream.
   *
   * The hydrodynamic fields, corresponding to density, velocity and stress, are
   * stored in LB_FluidNodes in the array lbfields, the populations in lbfluid
   * which is constructed as 2 x (Nx x Ny x Nz) x 19 array.
   */

/** Description of the LB Model in terms of the unit vectors of the 
 *  velocity sub-lattice and the corresponding coefficients 
 *  of the pseudo-equilibrium distribution */
typedef struct {

  /** number of velocities */
  int n_veloc ;

  /** unit vectors of the velocity sublattice */
  double (*c)[3];

  /** coefficients in the pseudo-equilibrium distribution */
  double (*coeff)[4];

  /** weights in the functional for the equilibrium distribution */
  double (*w);

  /** basis of moment space */
  double **e;

  /** speed of sound squared */
  double c_sound_sq;

} LB_Model;

/** Data structure for fluid on a local lattice site */
typedef struct {

  /** flag indicating whether fields have to be recomputed */
  int recalc_fields;

  /** local density */
  double rho[1];

  /** local momentum */
  double j[3];

  /** local stress tensor */
  double pi[6];

  /* local populations of the velocity directions
   *  are stored seperately to achieve higher performance */

  /** flag indicating whether a force is acting on the node */
  int has_force;

  /** local force density */
  double force[3];

#ifdef LB_BOUNDARIES
   /** flag indicating whether this site belongs to a boundary */
   int boundary;

  /** normal vector of the boundary surface */
  double *nvec; //doesn't work like that any more, I think (georg, 17.08.10)
#endif

} LB_FluidNode;

/** Data structure holding the parameters for the Lattice Boltzmann system. */
typedef struct {

  /** number density (LJ units) */
  double rho[LB_COMPONENTS];

  /** kinematic viscosity (LJ units) */
  double viscosity[LB_COMPONENTS];

  /** bulk viscosity (LJ units) */
  double bulk_viscosity[LB_COMPONENTS];

  /** lattice spacing (LJ units) */
  double agrid;

  /** time step for fluid propagation (LJ units)
   *  Note: Has to be larger than MD time step! */
  double tau;

  /** friction coefficient for viscous coupling (LJ units)
   * Note that the friction coefficient is quite high and may
   * lead to numerical artifacts with low order integrators */
  double friction[LB_COMPONENTS];

  /** external force applied to the fluid at each lattice site (MD units) */
  double ext_force[3]; /* Open question: Do we want a local force or global force? */
  double rho_lb_units[LB_COMPONENTS];
  double gamma_odd[LB_COMPONENTS];
  double gamma_even[LB_COMPONENTS];

  int resend_halo;
          
} LB_Parameters;

/** The DnQm model to be used. */
extern LB_Model lbmodel;

/** Struct holding the Lattice Boltzmann parameters */
extern LB_Parameters lbpar; 

/** The underlying lattice */
extern Lattice lblattice;

/** Pointer to the velocity populations of the fluid.
 * lbfluid[0] contains pre-collision populations, lbfluid[1]
 * contains post-collision populations*/
extern double **lbfluid[2];

/** Pointer to the hydrodynamic fields of the fluid */
extern LB_FluidNode *lbfields;

/** Switch indicating momentum exchange between particles and fluid */
extern int transfer_momentum;

/** Eigenvalue of collision operator corresponding to shear viscosity. */
extern double lblambda;

/** Eigenvalue of collision operator corresponding to bulk viscosity. */
extern double lblambda_bulk;

extern int resend_halo;

extern double gamma_shear;
extern double gamma_bulk;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Updates the Lattice Boltzmann system for one time step.
 * This function performs the collision step and the streaming step.
 * If external forces are present, they are applied prior to the collisions.
 * If boundaries are present, it also applies the boundary conditions.
 */
void lattice_boltzmann_update();

/** (Pre-)initializes data structures. */
void lb_pre_init();

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

/** Resets the forces on the fluid nodes */
void lb_reinit_forces();

/** Checks if all LB parameters are meaningful */
int lb_sanity_checks();

/** Sets the density and momentum on a local lattice site.
 * @param node  Pointer to the Node of the lattice site within the local domain (Input)
 * @param rho   Local density of the fluid (Input)
 * @param v     Local momentum of the fluid (Input)
 * @param pi    Local pressure of the fluid (Input)
 */
void lb_set_local_fields(LB_FluidNode *node, const double rho, const double *v, const double *pi);

/** Returns the mass, momentum and stress of a local lattice site.
 * @param node  The index of the lattice site within the local domain (Input)
 * @param rho   Local density of the fluid (Output)
 * @param j     Local momentum of the fluid (Output)
 * @param pi    Local stress tensor of the fluid (Output)
 */
void lb_get_local_fields(LB_FluidNode *node, double *rho, double *j, double *pi);

/** Calculates the equilibrium distributions.
    @param index Index of the local site
    @param rho local fluid density
    @param j local fluid speed
    @param pi local fluid pressure
*/
void lb_calc_n_from_rho_j_pi(const index_t index, const double rho, const double *j, double *pi);

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

/** calculates the fluid velocity at a given position of the 
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. */
int lb_lbfluid_get_interpolated_velocity(double* p, double* v); 

inline void lb_calc_local_fields(index_t index, double *rho, double *j, double *pi); 


/** Calculation of hydrodynamic modes.
 *
 *  @param index number of the node to calculate the modes for
 *  @param mode output pointer to a double[19] 
 */
void lb_calc_modes(index_t index, double *mode);

/** Calculate the local fluid density.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index the local lattice site (Input).
 * @param rho   local fluid density
 */
inline void lb_calc_local_rho(index_t index, double *rho) {

#ifndef D3Q19
#error Only D3Q19 is implemened!
#endif

  // unit conversion: mass density
  if (!(lattice_switch & LATTICE_LB)) {
    ERROR_SPRINTF(runtime_error(128), 
        "{ Error in lb_calc_local_rho in %s %d: CPU LB not switched on. } ", __FILE__, __LINE__);
    *rho =0;
    return;
  }

  double avg_rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid;

  *rho =   avg_rho
         + lbfluid[0][0][index]
         + lbfluid[0][1][index]  + lbfluid[0][2][index]
         + lbfluid[0][3][index]  + lbfluid[0][4][index]
         + lbfluid[0][5][index]  + lbfluid[0][6][index] 
         + lbfluid[0][7][index]  + lbfluid[0][8][index]  
	       + lbfluid[0][9][index]  + lbfluid[0][10][index]
         + lbfluid[0][11][index] + lbfluid[0][12][index] 
	       + lbfluid[0][13][index] + lbfluid[0][14][index] 
         + lbfluid[0][15][index] + lbfluid[0][16][index] 
	       + lbfluid[0][17][index] + lbfluid[0][18][index];

}

/** Calculate the local fluid momentum.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index The local lattice site (Input).
 * @param j local fluid speed
 */
inline void lb_calc_local_j(index_t index, double *j) {

#ifndef D3Q19
#error Only D3Q19 is implemened!
#endif
  if (!(lattice_switch & LATTICE_LB)) {
    ERROR_SPRINTF(runtime_error(128), 
        "{ Error in lb_calc_local_j in %s %d: CPU LB not switched on. } ", __FILE__, __LINE__);
    j[0]=j[1]=j[2]=0;
    return;
  }

  j[0] =   lbfluid[0][1][index]  - lbfluid[0][2][index]
         + lbfluid[0][7][index]  - lbfluid[0][8][index]  
         + lbfluid[0][9][index]  - lbfluid[0][10][index] 
         + lbfluid[0][11][index] - lbfluid[0][12][index] 
         + lbfluid[0][13][index] - lbfluid[0][14][index];
  j[1] =   lbfluid[0][3][index]  - lbfluid[0][4][index]
         + lbfluid[0][7][index]  - lbfluid[0][8][index]  
         - lbfluid[0][9][index]  + lbfluid[0][10][index]
         + lbfluid[0][15][index] - lbfluid[0][16][index] 
         + lbfluid[0][17][index] - lbfluid[0][18][index]; 
  j[2] =   lbfluid[0][5][index]  - lbfluid[0][6][index]  
         + lbfluid[0][11][index] - lbfluid[0][12][index] 
         - lbfluid[0][13][index] + lbfluid[0][14][index]
         + lbfluid[0][15][index] - lbfluid[0][16][index] 
         - lbfluid[0][17][index] + lbfluid[0][18][index];

}

/** Calculate the local fluid stress.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index The local lattice site (Input).
 * @param pi local fluid pressure
 */
inline void lb_calc_local_pi(index_t index, double *pi) {

  double rho;
  double j[3];
  
  if (!(lattice_switch & LATTICE_LB)) {
    ERROR_SPRINTF(runtime_error(128), 
        "{ Error in lb_calc_local_pi in %s %d: CPU LB not switched on. } ", __FILE__, __LINE__);
    j[0] = j[1] = j[2] = 0;
    return;
  }
  
  lb_calc_local_fields(index, &rho, j, pi);
}

/** Calculate the local fluid fields.
 * The calculation is implemented explicitly for the special case of D3Q19.
 *
 * @param index   Index of the local lattice site.
 * @param rho     local fluid density
 * @param j       local fluid speed
 * @param pi      local fluid pressure
 */
inline void lb_calc_local_fields(index_t index, double *rho, double *j, double *pi) {

  if (!(lattice_switch & LATTICE_LB)) {
    ERROR_SPRINTF(runtime_error(128), 
        "{ Error in lb_calc_local_fields in %s %d: CPU LB not switched on. } ", __FILE__, __LINE__);
    *rho=0; j[0]=j[1]=j[2]=0; pi[0]=pi[1]=pi[2]=pi[3]=pi[4]=pi[5]=0;
    return;
  }

#ifndef D3Q19
#error Only D3Q19 is implemened!
#endif
  
  if (!(lattice_switch & LATTICE_LB)) {
    ERROR_SPRINTF(runtime_error(128), 
        "{ Error in lb_calc_local_pi in %s %d: CPU LB not switched on. } ", __FILE__, __LINE__);
    j[0] = j[1] = j[2] = 0;
    return;
  }

#ifdef LB_BOUNDARIES
  if ( lbfields[index].boundary ) {
    *rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid;
    j[0] = 0.; j[1] = 0.;  j[2] = 0.;
    if (pi) {pi[0] = 0.; pi[1] = 0.; pi[2] = 0.; pi[3] = 0.; pi[4] = 0.; pi[5] = 0.;}
    return;
  }
#endif
  double mode[19];
  double modes_from_pi_eq[6];
  lb_calc_modes(index, mode);

  *rho = mode[0] + lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid;

  j[0] = mode[1];
  j[1] = mode[2];
  j[2] = mode[3];

#ifndef EXTERNAL_FORCES
  if (lbfields[index].has_force) 
#endif
  {
    j[0] += 0.5*lbfields[index].force[0];
    j[1] += 0.5*lbfields[index].force[1];
    j[2] += 0.5*lbfields[index].force[2];
  }
  if (!pi)
    return;

  /* equilibrium part of the stress modes */
  modes_from_pi_eq[0] = scalar(j,j)/ *rho;
  modes_from_pi_eq[1] = (SQR(j[0])-SQR(j[1]))/ *rho;
  modes_from_pi_eq[2] = (scalar(j,j) - 3.0*SQR(j[2]))/ *rho;
  modes_from_pi_eq[3] = j[0]*j[1]/ *rho;
  modes_from_pi_eq[4] = j[0]*j[2]/ *rho;
  modes_from_pi_eq[5] = j[1]*j[2]/ *rho;
  
  /* Now we must predict the outcome of the next collision */
  /* We immediately average pre- and post-collision. */
  mode[4] = modes_from_pi_eq[0] + (0.5+0.5*gamma_bulk )*(mode[4] - modes_from_pi_eq[0]);
  mode[5] = modes_from_pi_eq[1] + (0.5+0.5*gamma_shear)*(mode[5] - modes_from_pi_eq[1]);
  mode[6] = modes_from_pi_eq[2] + (0.5+0.5*gamma_shear)*(mode[6] - modes_from_pi_eq[2]);
  mode[7] = modes_from_pi_eq[3] + (0.5+0.5*gamma_shear)*(mode[7] - modes_from_pi_eq[3]);
  mode[8] = modes_from_pi_eq[4] + (0.5+0.5*gamma_shear)*(mode[8] - modes_from_pi_eq[4]);
  mode[9] = modes_from_pi_eq[5] + (0.5+0.5*gamma_shear)*(mode[9] - modes_from_pi_eq[5]);

  // Transform the stress tensor components according to the modes that
  // correspond to those used by U. Schiller. In terms of populations this
  // expression then corresponds exactly to those in Eqs. 116 - 121 in the
  // Duenweg and Ladd paper, when these are written out in populations.
  // But to ensure this, the expression in Schiller's modes has to be different!

  pi[0] = ( 2.0*(mode[0] + mode[4]) + mode[6] + 3.0*mode[5] )/6.0;  // xx
  pi[1] = mode[7];                                                  // xy
  pi[2] = ( 2.0*(mode[0] + mode[4]) + mode[6] - 3.0*mode[5] )/6.0;  // yy
  pi[3] = mode[8];                                                  // xz  
  pi[4] = mode[9];                                                  // yz
  pi[5] = ( mode[0] + mode[4] - mode[6] )/3.0;                      // zz

}

#ifdef LB_BOUNDARIES
inline void lb_local_fields_get_boundary_flag(index_t index, int *boundary) {
  
  if (!(lattice_switch & LATTICE_LB)) {
    ERROR_SPRINTF(runtime_error(128), 
        "{ Error in lb_local_fields_get_boundary_flag in %s %d: CPU LB not switched on. } ", __FILE__, __LINE__);
    *boundary = 0;
    return;
  }

  *boundary = lbfields[index].boundary;
}
#endif

/** Calculate the local fluid momentum.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index The local lattice site (Input).
 * @param pop fluid population
 */
inline void lb_get_populations(index_t index, double* pop) {
  int i=0;
  for (i=0; i<19; i++) {
    pop[i]=lbfluid[0][i][index]+lbmodel.coeff[i][0]*lbpar.rho[0];
  }
}

inline void lb_set_populations(index_t index, double* pop) {
  int i=0;
  for (i=0; i<19; i++) {
    lbfluid[0][i][index]=pop[i]-lbmodel.coeff[i][0]*lbpar.rho[0];
  }
}
#endif

#include "lbgpu.hpp"

#if defined (LB) || defined (LB_GPU)
/* A C level interface to the LB fluid */ 
int lb_lbfluid_set_density(double * p_dens);
int lb_lbfluid_set_visc(double * p_visc);
int lb_lbfluid_set_bulk_visc(double * p_bulk_visc);
int lb_lbfluid_set_gamma_odd(double * p_gamma_odd);
int lb_lbfluid_set_gamma_even(double * p_gamma_even);
int lb_lbfluid_set_friction(double * p_friction);
int lb_lbfluid_set_couple_flag(int couple_flag);
int lb_lbfluid_set_agrid(double p_agrid);
int lb_lbfluid_set_ext_force(int component, double p_fx, double p_fy, double p_fz);
int lb_lbfluid_set_tau(double p_tau);
int lb_lbfluid_set_remove_momentum(void);
#ifdef SHANCHEN
int lb_lbfluid_set_shanchen_coupling(double * p_coupling);
int lb_lbfluid_set_mobility(double * p_mobility);
#endif 

/* IO routines */
int lb_lbfluid_print_vtk_boundary(char* filename);
int lb_lbfluid_print_vtk_velocity(char* filename);
int lb_lbfluid_print_vtk_density(char** filename);
int lb_lbfluid_print_boundary(char* filename);
int lb_lbfluid_print_velocity(char* filename);

int lb_lbfluid_save_checkpoint(char* filename, int binary); 
int lb_lbfluid_load_checkpoint(char* filename, int binary);

int lb_lbnode_get_rho(int* ind, double* p_rho);
int lb_lbnode_get_u(int* ind, double* u);
int lb_lbnode_get_pi(int* ind, double* pi);
int lb_lbnode_get_pi_neq(int* ind, double* pi_neq);
int lb_lbnode_get_boundary(int* ind, int* p_boundary);
int lb_lbnode_get_pop(int* ind, double* pop);

int lb_lbnode_set_rho(int* ind, double *rho);
int lb_lbnode_set_u(int* ind, double* u);
int lb_lbnode_set_pi(int* ind, double* pi);
int lb_lbnode_set_pi_neq(int* ind, double* pi_neq);
int lb_lbnode_set_pop(int* ind, double* pop);

/** calculates the fluid velocity at a given position of the 
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. This version of the function
 * can be called without the position needing to be on the local processor */
int lb_lbfluid_get_interpolated_velocity_global(double* p, double* v); 

#endif

#endif /* _LB_H */

/*@}*/
