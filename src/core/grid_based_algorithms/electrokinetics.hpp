/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ELECTROKINETICS_HPP
#define _ELECTROKINETICS_HPP

#include "config.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"

// note that we need to declare the ek_parameters struct and instantiate it for
// LB_GPU to compile when electrokinetics is not compiled in. This seemed more
// elegant than ifdeffing multiple versions of the kernel integrate.
#ifdef CUDA

#if defined(EK_DOUBLE_PREC)
typedef double ekfloat;
#else
typedef float ekfloat;
#endif

#define MAX_NUMBER_OF_SPECIES 10

/* Data structure holding parameters and memory pointers for the link flux
 * system. */
typedef struct {
  float agrid;
  float time_step; // MD time step
  float lb_density;
  unsigned int dim_x;
  unsigned int dim_x_padded;
  unsigned int dim_y;
  unsigned int dim_z;
  unsigned int number_of_nodes;
  float viscosity;
  float bulk_viscosity;
  float gamma_odd;
  float gamma_even;
  float friction;
  float T;
  float prefactor;
  float lb_force_density[3];
  unsigned int number_of_species;
  int reaction_species[3];
  float rho_reactant_reservoir;
  float rho_product0_reservoir;
  float rho_product1_reservoir;
  float reaction_ct_rate;
  float reaction_fraction_0;
  float reaction_fraction_1;
  float mass_reactant;
  float mass_product0;
  float mass_product1;
  int stencil;
  int number_of_boundary_nodes;
  float fluctuation_amplitude;
  bool fluctuations;
  bool advection;
  bool fluidcoupling_ideal_contribution;
  bool es_coupling;
  float *charge_potential_buffer;
  float *electric_field;
  float *charge_potential;
  ekfloat *j;
  float *lb_force_density_previous;
#ifdef EK_DEBUG
  ekfloat *j_fluc;
#endif
  ekfloat *rho[MAX_NUMBER_OF_SPECIES];
  int species_index[MAX_NUMBER_OF_SPECIES];
  float density[MAX_NUMBER_OF_SPECIES];
  float D[MAX_NUMBER_OF_SPECIES];
  float d[MAX_NUMBER_OF_SPECIES];
  float valency[MAX_NUMBER_OF_SPECIES];
  float ext_force_density[3][MAX_NUMBER_OF_SPECIES];
  char *node_is_catalyst;
} EK_parameters;

#endif

#ifdef ELECTROKINETICS

/* Constants enumerating the links of a node in the link flux system EK_LINK_xyz
   is the number of the link in direction (x, y, z), where x, y and z can be 0,
   U or D representing 0 and one agrid in direction of or against the x, y or z
   axis. The numbering differs from the one used in the LB since the LB
   velocities are directed but the links are not. Links 0 - 8 represent
   the odd LB velocities and numbers 13 - 21 represent the even LB velocities
   (without the 0). In between there are the links connecting the corners, which
   represent the 3rd shell not used in the LB but in the advection. The
   following 13 constants are only defined for the sake of completeness.*/

#define EK_LINK_U00 0
#define EK_LINK_0U0 1
#define EK_LINK_00U 2
#define EK_LINK_UU0 3
#define EK_LINK_UD0 4
#define EK_LINK_U0U 5
#define EK_LINK_U0D 6
#define EK_LINK_0UU 7
#define EK_LINK_0UD 8

#define EK_LINK_UUU 9
#define EK_LINK_UUD 10
#define EK_LINK_UDU 11
#define EK_LINK_UDD 12

#define EK_LINK_D00 13
#define EK_LINK_0D0 14
#define EK_LINK_00D 15
#define EK_LINK_DD0 16
#define EK_LINK_DU0 17
#define EK_LINK_D0D 18
#define EK_LINK_D0U 19
#define EK_LINK_0DD 20
#define EK_LINK_0DU 21

#define EK_LINK_DDD 22
#define EK_LINK_DDU 23
#define EK_LINK_DUD 24
#define EK_LINK_DUU 25

extern EK_parameters ek_parameters;
extern bool ek_initialized;

void ek_integrate();
void ek_integrate_electrostatics();
void ek_print_parameters();
void ek_print_lbpar();
void lb_set_ek_pointer(EK_parameters *pointeradress);
unsigned int ek_calculate_boundary_mass();
int ek_print_vtk_density(int species, char *filename);
int ek_print_vtk_flux(int species, char *filename);
int ek_print_vtk_flux_fluc(int species, char *filename);
int ek_print_vtk_flux_link(int species, char *filename);
int ek_print_vtk_potential(char *filename);
int ek_print_vtk_particle_potential(char *filename);
int ek_print_vtk_lbforce_density(char *filename);
int ek_lb_print_vtk_density(char *filename);
int ek_lb_print_vtk_velocity(char *filename);
int ek_init();
int ek_set_agrid(float agrid);
int ek_set_lb_density(float lb_density);
int ek_set_viscosity(float viscosity);
int ek_set_friction(float friction);
int ek_set_T(float T);
int ek_set_prefactor(float prefactor);
int ek_set_electrostatics_coupling(bool electrostatics_coupling);
void ek_calculate_electrostatic_coupling();
int ek_set_bulk_viscosity(float bulk_viscosity);
int ek_set_gamma_odd(float gamma_odd);
int ek_set_gamma_even(float gamma_even);
int ek_set_density(int species, float density);
int ek_set_D(int species, float D);
int ek_set_valency(int species, float valency);
int ek_set_ext_force_density(int species, float ext_force_density_x,
                             float ext_force_density_y,
                             float ext_force_density_z);
int ek_set_stencil(int stencil);
int ek_set_advection(bool advection);
int ek_set_fluidcoupling(bool ideal_contribution);
int ek_set_fluctuations(bool fluctuations);
int ek_set_fluctuation_amplitude(float fluctuation_amplitude);
void ek_set_rng_state(uint64_t counter);
int ek_node_print_velocity(int x, int y, int z, double *velocity);
int ek_node_print_density(int species, int x, int y, int z, double *density);
int ek_node_print_flux(int species, int x, int y, int z, double *flux);
int ek_node_print_potential(int x, int y, int z, double *potential);
int ek_node_set_density(int species, int x, int y, int z, double density);
ekfloat ek_calculate_net_charge();
int ek_neutralize_system(int species);
int ek_save_checkpoint(char *filename, char *lb_filename);
int ek_load_checkpoint(char *filename);

#ifdef EK_BOUNDARIES
void ek_gather_wallcharge_species_density(ekfloat *wallcharge_species_density,
                                          int wallcharge_species);
void ek_init_species_density_wallcharge(ekfloat *wallcharge_species_density,
                                        int wallcharge_species);
#endif

#endif /* CUDA */

#endif /* ELECTROKINETICS_H */
