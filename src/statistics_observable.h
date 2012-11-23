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
#ifndef _STATISTICS_OBSERVABLE_H
#define _STATISTICS_OBSERVABLE_H

#include "config.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct s_observable;

struct s_observable {
  char* obs_name;
  void* container;
  int n;
  int (*update)    ( struct s_observable* obs );
  int (*calculate) ( struct s_observable* obs );
  double* last_value;
  int last_update;
};

typedef struct s_observable observable;

extern observable** observables;
extern int n_observables; 

int observable_calculate(observable* self);
int observable_update(observable* self);

/* Here we have the particular observables listed */
int observable_calc_particle_velocities(observable* self_);
int observable_calc_com_velocity(observable* self); 
int observable_calc_blocked_com_velocity(observable* self); 
/** Obtain the particle positions.
 * TODO: Folded or unfolded?
 */ 
int observable_calc_particle_positions(observable* self);
int observable_calc_particle_forces(observable* self);
int observable_calc_com_force(observable* self);
int observable_calc_blocked_com_force(observable* self);
int observable_calc_stress_tensor(observable* self);
int observable_calc_stress_tensor_acf_obs(observable* self);
int observable_calc_com_position(observable* self);
int observable_calc_blocked_com_position(observable* self);

#ifdef ELECTROSTATICS
int observable_calc_particle_currents(observable* self);
int observable_calc_currents(observable* self);
int observable_calc_dipole_moment(observable* self);
#endif

int observable_update_average(observable* self);
int observable_reset_average(observable* self);
typedef struct {
  observable* reference_observable;
  unsigned int n_sweeps;
} observable_average_container;

/** Calculate structure factor from positions and scattering length */
int observable_calc_structure_factor(observable* self);
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
int observable_calc_interacts_with(observable* self);
typedef struct {
  double cutoff;
  IntList *ids1;
  IntList *ids2;
} iw_params;


/** Do nothing */
int observable_calc_obs_nothing (observable* self);

int observable_calc_flux_density_profile(observable* self);
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

int observable_calc_density_profile(observable* self);

int observable_calc_lb_velocity_profile(observable* self);

int observable_calc_radial_density_profile(observable* self);
int observable_calc_radial_flux_density_profile(observable* self);
int observable_calc_lb_radial_velocity_profile(observable* self);
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
