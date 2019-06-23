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
#include "lb-d3q19.hpp"

#ifdef CUDA

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
LB_rho_v_pi_gpu *host_values = nullptr;

static int max_ran = 1000000;
// static double tau;

/** measures the MD time since the last fluid update */
static int fluidstep = 0;

// clock_t start, end;
int i;

int n_extern_node_force_densities = 0;
LB_extern_nodeforcedensity_gpu *host_extern_node_force_densities = nullptr;
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

/** Release the fluid. */
/*not needed in Espresso but still not deleted.
  Despite the name (TODO: change it), it releases
  only the fluid-related memory on the gpu.*/
void lb_release_gpu() {

  if (host_values != nullptr) {
    free(host_values);
    host_values = nullptr;
  }
  //  if(host_forces!=nullptr) free(host_forces);
  //  if(host_data  !=nullptr) free(host_data);
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
        box_l[0] / lbpar_gpu.agrid); // TODO code duplication with lb.c start
    lbpar_gpu.dim_y = (unsigned int)round(box_l[1] / lbpar_gpu.agrid);
    lbpar_gpu.dim_z = (unsigned int)round(box_l[2] / lbpar_gpu.agrid);

    unsigned int tmp[3];

    tmp[0] = lbpar_gpu.dim_x;
    tmp[1] = lbpar_gpu.dim_y;
    tmp[2] = lbpar_gpu.dim_z;

    /* sanity checks */
    int dir;

    for (dir = 0; dir < 3; dir++) {
      /* check if box_l is compatible with lattice spacing */
      if (fabs(box_l[dir] - tmp[dir] * lbpar_gpu.agrid) > 1.0e-3) {
        runtimeErrorMsg() << "Lattice spacing lbpar_gpu.agrid= "
                          << lbpar_gpu.agrid << " is incompatible with box_l["
                          << dir << "]=" << box_l[dir];
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

int lb_lbnode_set_extforce_density_GPU(int const ind[3], double const f[3]) {
  if (ind[0] < 0 || ind[0] >= int(lbpar_gpu.dim_x) || ind[1] < 0 ||
      ind[1] >= int(lbpar_gpu.dim_y) || ind[2] < 0 ||
      ind[2] >= int(lbpar_gpu.dim_z))
    return ES_ERROR;

  unsigned int index = ind[0] + ind[1] * lbpar_gpu.dim_x +
                       ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;

  size_t size_of_extforces = (n_extern_node_force_densities + 1) *
                             sizeof(LB_extern_nodeforcedensity_gpu);
  host_extern_node_force_densities =
      (LB_extern_nodeforcedensity_gpu *)Utils::realloc(
          host_extern_node_force_densities, size_of_extforces);

  host_extern_node_force_densities[n_extern_node_force_densities]
      .force_density[0] = (float)f[0];
  host_extern_node_force_densities[n_extern_node_force_densities]
      .force_density[1] = (float)f[1];
  host_extern_node_force_densities[n_extern_node_force_densities]
      .force_density[2] = (float)f[2];

  host_extern_node_force_densities[n_extern_node_force_densities].index = index;
  n_extern_node_force_densities++;

  if (lbpar_gpu.external_force_density == 0)
    lbpar_gpu.external_force_density = 1;

  lb_init_extern_nodeforcedensities_GPU(n_extern_node_force_densities,
                                        host_extern_node_force_densities,
                                        &lbpar_gpu);

  return ES_OK;
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

void lb_lbfluid_particles_add_momentum(float const momentum[3]) {
  auto &parts = partCfg();
  auto const n_part = parts.size();

  // set_particle_v invalidates the parts pointer, so we need to defer setting
  // the new values
  std::vector<std::pair<int, double[3]>> new_velocity(n_part);

  size_t i = 0;
  for (auto const &p : parts) {
    new_velocity[i].first = p.p.identity;
    const auto factor = 1 / (p.p.mass * n_part);
    new_velocity[i].second[0] = p.m.v[0] + momentum[0] * factor;
    new_velocity[i].second[1] = p.m.v[1] + momentum[1] * factor;
    new_velocity[i].second[2] = p.m.v[2] + momentum[2] * factor;
    ++i;
  }
  for (auto &p : new_velocity) {
    set_particle_v(p.first, p.second);
  }
}

void lb_lbfluid_calc_linear_momentum(float momentum[3], int include_particles,
                                     int include_lbfluid) {
  auto linear_momentum =
      calc_linear_momentum(include_particles, include_lbfluid);
  momentum[0] = linear_momentum[0];
  momentum[1] = linear_momentum[1];
  momentum[2] = linear_momentum[2];
}

void lb_set_agrid_gpu(double agrid) {
  lbpar_gpu.agrid = static_cast<float>(agrid);

  lbpar_gpu.dim_x = static_cast<unsigned int>(rint(box_l[0] / agrid));
  lbpar_gpu.dim_y = static_cast<unsigned int>(rint(box_l[1] / agrid));
  lbpar_gpu.dim_z = static_cast<unsigned int>(rint(box_l[2] / agrid));
  unsigned int tmp[3];
  tmp[0] = lbpar_gpu.dim_x;
  tmp[1] = lbpar_gpu.dim_y;
  tmp[2] = lbpar_gpu.dim_z;
  /* sanity checks */
  for (int dir = 0; dir < 3; dir++) {
    /* check if box_l is compatible with lattice spacing */
    if (fabs(box_l[dir] - tmp[dir] * agrid) > ROUND_ERROR_PREC) {
      runtimeErrorMsg() << "Lattice spacing p_agrid= " << agrid
                        << " is incompatible with box_l[" << dir
                        << "]=" << box_l[dir] << ", factor=" << tmp[dir]
                        << " err= " << fabs(box_l[dir] - tmp[dir] * agrid);
    }
  }
  lbpar_gpu.number_of_nodes =
      lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
}

#endif /*  CUDA */
