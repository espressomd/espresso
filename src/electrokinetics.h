/*
   Copyright (C) 2010,2011,2012 The ESPResSo project

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

#ifndef ELECTROKINETICS_H
#define ELECTROKINETICS_H

#include <cufft.h>
#include "config.h"


#ifdef __cplusplus
extern "C" {
#endif

//#ifdef ELECTROKINETICS

#define MAX_NUMBER_OF_SPECIES 10

//TODO update

/* Constants enumerating the links of a node in the link flux system EK_LINK_xyz
   is the number of the link in direction (x, y, z), where x, y and z can be 0, 
   U or D representing 0 and one agrid in direction of or against the x, y or z 
   axis. The numbering differs from the one used in the LB since the LB 
   velocities are directed but the links are not. The first 9 links represent 
   the odd LB velocities and the last nine represent the even LB velocities 
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

/* Data structure holding parameters and memory pointers for the link flux system. */

typedef struct {
  float agrid;
  float time_step; //MD time step
  unsigned int dim_x;
  unsigned int dim_y;
  unsigned int dim_z;
  unsigned int number_of_nodes;
  float viscosity;
  float bulk_viscosity;
  float gamma_odd;
  float gamma_even;
  float friction;
  float T;
  float bjerrumlength;
  unsigned int number_of_species;
  unsigned int accelerated_frame_enabled;
  float accelerated_frame_boundary_mass;
  int reaction_species[3];
  float rho_reactant_reservoir;
  float rho_product0_reservoir;
  float rho_product1_reservoir;
  float reaction_ct_rate;
  float reaction_radius;
  float reaction_fraction_0;
  float reaction_fraction_1;
  cufftReal* greensfcn;
  cufftComplex* charge_potential;
  float* j;
  float* rho[MAX_NUMBER_OF_SPECIES];
  int species_index[MAX_NUMBER_OF_SPECIES];
  float density[MAX_NUMBER_OF_SPECIES];
  float D[MAX_NUMBER_OF_SPECIES];
  float d[MAX_NUMBER_OF_SPECIES];
  float valency[MAX_NUMBER_OF_SPECIES];
  float ext_force[3][MAX_NUMBER_OF_SPECIES];
  char* node_is_catalyst;
#ifdef EK_REACTION
  float* pressure;
#endif
} EK_parameters;

extern EK_parameters ek_parameters;

void ek_integrate();
void ek_print_parameters();
void ek_print_lbpar();
int ek_print_vtk_density(int species, char* filename);
int ek_print_vtk_flux(int species, char* filename);
int ek_print_vtk_potential(char* filename);
int ek_print_vtk_lbforce(char* filename);
int ek_lb_print_vtk_density(char* filename);
int ek_lb_print_vtk_velocity(char* filename);
int ek_init();
int ek_set_agrid(double agrid);
int ek_set_viscosity(double viscosity);
int ek_set_friction(double friction);
int ek_set_T(double T);
int ek_set_bjerrumlength(double bjerrumlength);
int ek_set_bulk_viscosity(double bulk_viscosity);
int ek_set_gamma_odd(double gamma_odd);
int ek_set_gamma_even(double gamma_even);
int ek_set_density(int species, double density);
int ek_set_D(int species, double D);
int ek_set_valency(int species, double valency);
int ek_set_ext_force(int species, double ext_force_x, double ext_force_y, double ext_force_z);
int ek_set_accelerated_frame( int enabled, double boundary_mass );

#ifdef EK_BOUNDARIES
void ek_init_species_density_wallcharge(float* wallcharge_species_density, int wallcharge_species);
#endif


#ifdef EK_REACTION
int ek_set_reaction(int reactant, int product0, int product1, float rho_reactant_reservoir, float rho_product0_reservoir, float rho_product1_reservoir, float reaction_ct_rate, float reaction_radius, float reaction_fraction_0, float reaction_fraction_1 );
int ek_print_vtk_pressure(char* filename);
#endif

//#endif /* ELECTROKINETICS */

#ifdef __cplusplus
}
#endif

#endif /* ELECTROKINETICS_H */
