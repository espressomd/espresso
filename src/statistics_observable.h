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
#ifndef _STATISTICS_OBSERVABLE_H
#define _STATISTICS_OBSERVABLE_H

#include "config.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct {
  char* obs_name;
  void* args;
  int n;
  int (*fun)  ( void* A_args, double* A, unsigned int dim_A);
  double* last_value;
  int last_update;
} observable;

extern observable** observables;
extern int n_observables; 

int* observable_calc(observable* self, double* A);

/* Here we have the particular observables listed */
int observable_particle_velocities(void* idlist, double* A, unsigned int n_A);
int observable_particle_angular_momentum(void* idlist, double* A, unsigned int n_A);
int observable_com_velocity(void* idlist, double* A, unsigned int n_A); 
int observable_blocked_com_velocity(void* idlist, double* A, unsigned int n_A); 
/** Obtain the particle positions.
 * TODO: Folded or unfolded?
 */ 
int observable_particle_positions(void* typelist, double* A, unsigned int n_A);
int observable_particle_forces(void* typelist, double* A, unsigned int n_A);
int observable_com_force(void* typelist, double* A, unsigned int n_A);
int observable_blocked_com_force(void* typelist, double* A, unsigned int n_A);
int observable_stress_tensor(void* typelist, double* A, unsigned int n_A);
int observable_stress_tensor_acf_obs(void* typelist, double* A, unsigned int n_A);
int observable_com_position(void* idlist, double* A, unsigned int n_A);
int observable_blocked_com_position(void* idlist, double* A, unsigned int n_A);

#ifdef ELECTROSTATICS
int observable_particle_currents(void* typelist, double* A, unsigned int n_A);
int observable_currents(void* typelist, double* A, unsigned int n_A);
int observable_dipole_moment(void* typelist, double* A, unsigned int n_A);
#endif


/** Calculate structure factor from positions and scattering length */
int observable_structure_factor(void* params, double* A, unsigned int n_A);
typedef struct {
// FIXME finish the implementation of scattering length
  IntList* id_list;
  DoubleList *scattering_length; // Scattering lengths of particles
  int order;
  int dim_sf; // number of q vectors
  int *q_vals; // values of q vectors
  double *q_density; // number of q vectors per bin
  // entries for spherical averaging
} observable_sf_params;

/** See if particles from idList1 interact with any of the particles in idList2 
input parameters are passed via struct iw_params
*/
int observable_interacts_with(void* params, double* A, unsigned int n_A);
typedef struct {
  double cutoff;
  IntList *ids1;
  IntList *ids2;
} iw_params;


/** Do nothing */
int observable_obs_nothing (void* params, double* A, unsigned int n_A);

int observable_flux_density_profile(void* params, double* A, unsigned int n_A);
typedef struct { 
  IntList* id_list;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
  int xbins;
  int ybins;
  int zbins;
} profile_data;

int observable_density_profile(void* params, double* A, unsigned int n_A);

int observable_lb_velocity_profile(void* params, double* A, unsigned int n_A);

int observable_radial_density_profile(void* params, double* A, unsigned int n_A);
int observable_radial_flux_density_profile(void* params, double* A, unsigned int n_A);
int observable_lb_radial_velocity_profile(void* params, double* A, unsigned int n_A);
typedef struct {
  IntList* id_list;
  double minr;
  double maxr;
  double minphi;
  double maxphi;
  double minz;
  double maxz;
  double center[3];
  double axis[3];
  int phibins;
  int rbins;
  int zbins;
} radial_profile_data;



#endif
