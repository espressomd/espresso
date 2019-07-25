/*
   Copyright (C) 2010-2018 The ESPResSo project

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
/** \file
 *  %Lattice Boltzmann on GPUs.
 *
 *  The corresponding header file is lbgpu.hpp.
 */

#include "lbgpu.hpp"
#ifdef CUDA
#include "errorhandling.hpp"
#include "lb-d3q19.hpp"

#include "communication.hpp"
#include "cuda_interface.hpp"
#include "debug.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lbgpu.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "statistics.hpp"

#include <utils/constants.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <random>

LB_particle_allocation_state lb_reinit_particles_gpu;

LB_parameters_gpu lbpar_gpu = {
    // rho
    0.0,
    // mu
    0.0,
    // viscosity
    0.0,
    // gamma_shear
    0.0,
    // gamma_bulk
    0.0,
    // gamma_odd
    0.0,
    // gamma_even
    0.0,
    // is_TRT
    false,
    // bulk_viscosity
    -1.0,
    // agrid
    -1.0,
    // tau
    -1.0,
    // time_step
    0.0,
    // dim_x;
    0,
    // dim_y;
    0,
    // dim_z;
    0,
    // number_of_nodes
    0,
    // number_of_particles
    0,
#ifdef LB_BOUNDARIES_GPU
    // number_of_boundnodes
    0,
#endif
    // calc_val
    1,
    // external_force
    0,
    // ext_force
    {0.0, 0.0, 0.0},
    // reinit
    0,
    // Thermal energy
    0.0};

/** this is the array that stores the hydrodynamic fields for the output */
std::vector<LB_rho_v_pi_gpu> host_values(0);

static int max_ran = 1000000;
// static double tau;

/** measures the MD time since the last fluid update */
static int fluidstep = 0;

// clock_t start, end;
int i;

int n_extern_node_force_densities = 0;
std::vector<LB_extern_nodeforcedensity_gpu> host_extern_node_force_densities;
bool ek_initialized = false;

/*-----------------------------------------------------------*/
/** main of lb_gpu_programm */
/*-----------------------------------------------------------*/

/** %Lattice Boltzmann update gpu called from integrate.cpp */
void lattice_boltzmann_update_gpu() {

  auto factor = (int)round(lbpar_gpu.tau / time_step);

  fluidstep += 1;

  if (fluidstep >= factor) {

    fluidstep = 0;
    lb_integrate_GPU();
    LB_TRACE(fprintf(stderr, "lb_integrate_GPU \n"));
  }
}

/** (Re-)allocation of the memory needed for the particles (CPU part) */
void lb_realloc_particles_gpu() {

  lbpar_gpu.number_of_particles = n_part;
  LB_TRACE(printf("#particles realloc\t %u \n", lbpar_gpu.number_of_particles));

  lb_realloc_particles_GPU_leftovers(&lbpar_gpu);
}

/** (Re-)initialize the fluid according to the given value of rho. */
void lb_reinit_fluid_gpu() {

  lb_reinit_parameters_gpu();
  if (lbpar_gpu.number_of_nodes != 0) {
    lb_reinit_GPU(&lbpar_gpu);
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
    lbpar_gpu.reinit = 1;
  }

  LB_TRACE(fprintf(stderr, "lb_reinit_fluid_gpu \n"));
}

/** (Re-)initialize the fluid. */
void lb_reinit_parameters_gpu() {
  lbpar_gpu.time_step = (float)time_step;
  lbpar_gpu.mu = 0.0;

  if (lbpar_gpu.viscosity > 0.0 && lbpar_gpu.agrid > 0.0 &&
      lbpar_gpu.tau > 0.0) {
    /* Eq. (80) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar_gpu.gamma_shear = 1. - 2. / (6. * lbpar_gpu.viscosity + 1.);
  }

  if (lbpar_gpu.bulk_viscosity > 0.0) {
    /* Eq. (81) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar_gpu.gamma_bulk = 1. - 2. / (9. * lbpar_gpu.bulk_viscosity + 1.);
  }

  // By default, gamma_even and gamma_odd are chosen such that the MRT becomes
  // a TRT with ghost mode relaxation factors that minimize unphysical wall
  // slip at bounce-back boundaries. For the relation between the gammas
  // achieving this, consult
  //  D. d’Humières, I. Ginzburg, Comp. & Math. w. App. 58(5):823–840 (2009)
  // Note that the relaxation operator in Espresso is defined as
  //  m* = m_eq + gamma * (m - m_eq)
  // as opposed to this reference, where
  //  m* = m + lambda * (m - m_eq)

  if (lbpar_gpu.is_TRT) {
    lbpar_gpu.gamma_bulk = lbpar_gpu.gamma_shear;
    lbpar_gpu.gamma_even = lbpar_gpu.gamma_shear;
    lbpar_gpu.gamma_odd =
        -(7.0f * lbpar_gpu.gamma_even + 1.0f) / (lbpar_gpu.gamma_even + 7.0f);
  }

  if (lbpar_gpu.kT > 0.0) { /* fluctuating hydrodynamics ? */

    LB_TRACE(fprintf(stderr, "fluct on \n"));
    /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).*/
    /* Note that the modes are not normalized as in the paper here! */
    lbpar_gpu.mu = lbpar_gpu.kT * lbpar_gpu.tau * lbpar_gpu.tau /
                   D3Q19::c_sound_sq<float> /
                   (lbpar_gpu.agrid * lbpar_gpu.agrid);
  }
  LB_TRACE(fprintf(stderr, "lb_reinit_prarameters_gpu \n"));

#ifdef ELECTROKINETICS
  if (ek_initialized) {
    lbpar_gpu.dim_x = (unsigned int)round(
        box_geo.length()[0] /
        lbpar_gpu.agrid); // TODO code duplication with lb.c start
    lbpar_gpu.dim_y =
        (unsigned int)round(box_geo.length()[1] / lbpar_gpu.agrid);
    lbpar_gpu.dim_z =
        (unsigned int)round(box_geo.length()[2] / lbpar_gpu.agrid);

    unsigned int tmp[3];

    tmp[0] = lbpar_gpu.dim_x;
    tmp[1] = lbpar_gpu.dim_y;
    tmp[2] = lbpar_gpu.dim_z;

    /* sanity checks */
    int dir;

    for (dir = 0; dir < 3; dir++) {
      /* check if box_l is compatible with lattice spacing */
      if (fabs(box_geo.length()[dir] - tmp[dir] * lbpar_gpu.agrid) > 1.0e-3) {
        runtimeErrorMsg() << "Lattice spacing lbpar_gpu.agrid= "
                          << lbpar_gpu.agrid << " is incompatible with box_l["
                          << dir << "]=" << box_geo.length()[dir];
      }
    }

    lbpar_gpu.number_of_nodes =
        lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
    lbpar_gpu.tau = (float)time_step; // TODO code duplication with lb.c end
  }
#endif

  LB_TRACE(fprintf(stderr, "lb_reinit_prarameters_gpu \n"));

  reinit_parameters_GPU(&lbpar_gpu);
}

/** Performs a full initialization of the lattice Boltzmann system.
 *  All derived parameters and the fluid are reset to their default values.
 */
void lb_init_gpu() {

  LB_TRACE(printf("Begin initializing fluid on GPU\n"));
  /** set parameters for transfer to gpu */
  lb_reinit_parameters_gpu();

  lb_realloc_particles_gpu();

  lb_init_GPU(&lbpar_gpu);

  gpu_init_particle_comm();
  cuda_bcast_global_part_params();

  LB_TRACE(printf("Initializing fluid on GPU successful\n"));
}

void lb_GPU_sanity_checks() {
  if (this_node == 0) {
    if (lbpar_gpu.agrid < 0.0) {
      runtimeErrorMsg() << "Lattice Boltzmann agrid not set";
    }
    if (lbpar_gpu.tau < 0.0) {
      runtimeErrorMsg() << "Lattice Boltzmann time step not set";
    }
    if (lbpar_gpu.rho < 0.0) {
      runtimeErrorMsg() << "Lattice Boltzmann fluid density not set";
    }
    if (lbpar_gpu.viscosity < 0.0) {
      runtimeErrorMsg() << "Lattice Boltzmann fluid viscosity not set";
    }
  }
}

void lb_set_agrid_gpu(double agrid) {
  lbpar_gpu.agrid = static_cast<float>(agrid);

  lbpar_gpu.dim_x =
      static_cast<unsigned int>(rint(box_geo.length()[0] / agrid));
  lbpar_gpu.dim_y =
      static_cast<unsigned int>(rint(box_geo.length()[1] / agrid));
  lbpar_gpu.dim_z =
      static_cast<unsigned int>(rint(box_geo.length()[2] / agrid));
  unsigned int tmp[3];
  tmp[0] = lbpar_gpu.dim_x;
  tmp[1] = lbpar_gpu.dim_y;
  tmp[2] = lbpar_gpu.dim_z;
  /* sanity checks */
  for (int dir = 0; dir < 3; dir++) {
    /* check if box_l is compatible with lattice spacing */
    if (fabs(box_geo.length()[dir] - tmp[dir] * agrid) > ROUND_ERROR_PREC) {
      runtimeErrorMsg() << "Lattice spacing p_agrid= " << agrid
                        << " is incompatible with box_l[" << dir
                        << "]=" << box_geo.length()[dir]
                        << ", factor=" << tmp[dir] << " err= "
                        << fabs(box_geo.length()[dir] - tmp[dir] * agrid);
    }
  }
  lbpar_gpu.number_of_nodes =
      lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
}

#endif /*  CUDA */
