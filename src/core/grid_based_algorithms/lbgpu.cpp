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
/** \file
 *  %Lattice Boltzmann on GPUs.
 *
 *  The corresponding header file is lbgpu.hpp.
 */

#include "lbgpu.hpp"
#ifdef CUDA
#include "errorhandling.hpp"
#include "lb-d3q19.hpp"

#include "Particle.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lbgpu.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "statistics.hpp"

#include <utils/constants.hpp>

#include <cmath>
#include <cstdlib>

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

/** measures the MD time since the last fluid update */
static int fluidstep = 0;

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
  }
}

/** (Re-)initialize the fluid according to the given value of rho. */
void lb_reinit_fluid_gpu() {

  lb_reinit_parameters_gpu();
  if (lbpar_gpu.number_of_nodes != 0) {
    lb_reinit_GPU(&lbpar_gpu);
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
    lbpar_gpu.reinit = 1;
  }
}

/** (Re-)initialize the fluid.
 *  See @cite dunweg07a and @cite dhumieres09a.
 */
void lb_reinit_parameters_gpu() {
  lbpar_gpu.time_step = static_cast<float>(time_step);
  lbpar_gpu.mu = 0.0;

  if (lbpar_gpu.viscosity > 0.0 && lbpar_gpu.agrid > 0.0 &&
      lbpar_gpu.tau > 0.0) {
    /* Eq. (80) @cite dunweg07a. */
    lbpar_gpu.gamma_shear = 1.f - 2.f / (6.f * lbpar_gpu.viscosity + 1.f);
  }

  if (lbpar_gpu.bulk_viscosity > 0.0) {
    /* Eq. (81) @cite dunweg07a. */
    lbpar_gpu.gamma_bulk = 1.f - 2.f / (9.f * lbpar_gpu.bulk_viscosity + 1.f);
  }

  // By default, gamma_even and gamma_odd are chosen such that the MRT becomes
  // a TRT with ghost mode relaxation factors that minimize unphysical wall
  // slip at bounce-back boundaries. For the relation between the gammas
  // achieving this, consult @cite dhumieres09a.
  // Note that the relaxation operator in ESPResSo is defined as
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

    /* Eq. (51) @cite dunweg07a.*/
    /* Note that the modes are not normalized as in the paper here! */
    lbpar_gpu.mu = lbpar_gpu.kT * lbpar_gpu.tau * lbpar_gpu.tau /
                   D3Q19::c_sound_sq<float> /
                   (lbpar_gpu.agrid * lbpar_gpu.agrid);
  }

#ifdef ELECTROKINETICS
  if (ek_initialized) {
    lbpar_gpu.dim_x =
        static_cast<unsigned int>(round(box_geo.length()[0] / lbpar_gpu.agrid));
    lbpar_gpu.dim_y =
        static_cast<unsigned int>(round(box_geo.length()[1] / lbpar_gpu.agrid));
    lbpar_gpu.dim_z =
        static_cast<unsigned int>(round(box_geo.length()[2] / lbpar_gpu.agrid));

    unsigned int tmp[3];
    tmp[0] = lbpar_gpu.dim_x;
    tmp[1] = lbpar_gpu.dim_y;
    tmp[2] = lbpar_gpu.dim_z;

    /* sanity checks */
    for (int dir = 0; dir < 3; dir++) {
      /* check if box_l is compatible with lattice spacing */
      if (fabs(box_geo.length()[dir] -
               static_cast<float>(tmp[dir]) * lbpar_gpu.agrid) > 1.0e-3) {
        runtimeErrorMsg() << "Lattice spacing lbpar_gpu.agrid= "
                          << lbpar_gpu.agrid << " is incompatible with box_l["
                          << dir << "]=" << box_geo.length()[dir];
      }
    }

    lbpar_gpu.number_of_nodes =
        lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
    lbpar_gpu.tau = static_cast<float>(time_step);
  }
#endif

  reinit_parameters_GPU(&lbpar_gpu);
}

/** Performs a full initialization of the lattice Boltzmann system.
 *  All derived parameters and the fluid are reset to their default values.
 */
void lb_init_gpu() {
  /* set parameters for transfer to gpu */
  lb_reinit_parameters_gpu();

  lb_init_GPU(&lbpar_gpu);

  gpu_init_particle_comm();
  cuda_bcast_global_part_params();
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
      static_cast<unsigned int>(round(box_geo.length()[0] / agrid));
  lbpar_gpu.dim_y =
      static_cast<unsigned int>(round(box_geo.length()[1] / agrid));
  lbpar_gpu.dim_z =
      static_cast<unsigned int>(round(box_geo.length()[2] / agrid));
  unsigned int tmp[3];
  tmp[0] = lbpar_gpu.dim_x;
  tmp[1] = lbpar_gpu.dim_y;
  tmp[2] = lbpar_gpu.dim_z;
  /* sanity checks */
  for (int dir = 0; dir < 3; dir++) {
    /* check if box_l is compatible with lattice spacing */
    if (fabs(box_geo.length()[dir] - static_cast<float>(tmp[dir]) * agrid) >
        ROUND_ERROR_PREC) {
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
