/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
 *  %Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 *  Includes fluctuating LB and coupling to MD particles via frictional
 *  momentum transfer.
 *
 *  The corresponding header file is lb.hpp.
 */

#include "grid_based_algorithms/lb.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "debug.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "halo.hpp"
#include "integrate.hpp"
#include "lb-d3q19.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "random.hpp"
#include "virtual_sites/lb_inertialess_tracers.hpp"

#include "utils/u32_to_u64.hpp"
#include <utils/Counter.hpp>
#include <utils/index.hpp>
#include <utils/math/matrix_vector_product.hpp>
#include <utils/uniform.hpp>
using Utils::get_linear_index;
#include <utils/constants.hpp>

#include <Random123/philox.h>
#include <boost/multi_array.hpp>
#include <mpi.h>
#include <profiler/profiler.hpp>

#include <cassert>
#include <cinttypes>
#include <cstdio>
#include <fstream>
#include <iostream>

namespace {
/** Basis of the mode space as described in [Duenweg, Schiller, Ladd] */
extern constexpr const std::array<std::array<int, 19>, 19> e_ki = {
    {{{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
     {{0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0}},
     {{0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1}},
     {{0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}},
     {{-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
     {{0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1}},
     {{-0, 1, 1, 1, 1, -2, -2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1}},
     {{0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0}},
     {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0}},
     {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1}},
     {{0, -2, 2, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0}},
     {{0, 0, 0, -2, 2, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1}},
     {{0, 0, 0, 0, 0, -2, 2, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}},
     {{0, -0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0}},
     {{0, 0, 0, -0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, -1, 1, -1, 1}},
     {{0, 0, 0, 0, 0, -0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1}},
     {{1, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
     {{0, -1, -1, 1, 1, -0, -0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1}},
     {{0, -1, -1, -1, -1, 2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1}}}};

extern constexpr const std::array<std::array<int, 19>, 19> e_ki_transposed = {
    {{{1, 0, 0, 0, -1, 0, -0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}},
     {{1, 1, 0, 0, 0, 1, 1, 0, 0, 0, -2, 0, 0, -0, 0, 0, -2, -1, -1}},
     {{1, -1, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, -1, -1}},
     {{1, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0, -2, 0, 0, -0, 0, -2, 1, -1}},
     {{1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, -2, 1, -1}},
     {{1, 0, 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, -2, 0, 0, -0, -2, -0, 2}},
     {{1, 0, 0, -1, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, -0, 2}},
     {{1, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 2}},
     {{1, -1, -1, 0, 1, 0, 2, 1, 0, 0, -1, -1, 0, -1, -1, 0, 1, 0, 2}},
     {{1, 1, -1, 0, 1, 0, 2, -1, 0, 0, 1, -1, 0, 1, -1, 0, 1, 0, 2}},
     {{1, -1, 1, 0, 1, 0, 2, -1, 0, 0, -1, 1, 0, -1, 1, 0, 1, 0, 2}},
     {{1, 1, 0, 1, 1, 1, -1, 0, 1, 0, 1, 0, 1, -1, 0, 1, 1, 1, -1}},
     {{1, -1, 0, -1, 1, 1, -1, 0, 1, 0, -1, 0, -1, 1, 0, -1, 1, 1, -1}},
     {{1, 1, 0, -1, 1, 1, -1, 0, -1, 0, 1, 0, -1, -1, 0, -1, 1, 1, -1}},
     {{1, -1, 0, 1, 1, 1, -1, 0, -1, 0, -1, 0, 1, 1, 0, 1, 1, 1, -1}},
     {{1, 0, 1, 1, 1, -1, -1, 0, 0, 1, 0, 1, 1, 0, -1, -1, 1, -1, -1}},
     {{1, 0, -1, -1, 1, -1, -1, 0, 0, 1, 0, -1, -1, 0, 1, 1, 1, -1, -1}},
     {{1, 0, 1, -1, 1, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 1, -1, -1}},
     {{1, 0, -1, 1, 1, -1, -1, 0, 0, -1, 0, -1, 1, 0, 1, -1, 1, -1, -1}}}};
} // namespace

#ifdef ADDITIONAL_CHECKS
static void lb_check_halo_regions(const LB_Fluid &lbfluid);
#endif // ADDITIONAL_CHECKS

/** Counter for the RNG */
boost::optional<Utils::Counter<uint64_t>> rng_counter_fluid;

/** Struct holding the Lattice Boltzmann parameters */
LB_Parameters lbpar = {
    // density
    0.0,
    // viscosity
    0.0,
    // bulk_viscosity
    -1.0,
    // agrid
    -1.0,
    // tau
    -1.0,
    // ext_force_density
    {0.0, 0.0, 0.0},
    // gamma_odd
    0.,
    // gamma_even
    0.,
    // gamma_shear
    0.,
    // gamma_bulk
    0.,
    // is_TRT
    false,
    // phi
    {},
    // Thermal energy
    0.0};

/** The underlying lattice structure */
Lattice lblattice;

using LB_FluidData = boost::multi_array<double, 2>;
static LB_FluidData lbfluid_a;
static LB_FluidData lbfluid_b;

/** Pointer to the velocity populations of the fluid.
 *  lbfluid contains pre-collision populations, lbfluid_post
 *  contains post-collision.
 */
LB_Fluid lbfluid;
LB_Fluid lbfluid_post;

/** Pointer to the hydrodynamic fields of the fluid nodes */
std::vector<LB_FluidNode> lbfields;

/** Communicator for halo exchange between processors */
HaloCommunicator update_halo_comm = {0, nullptr};

/** measures the MD time since the last fluid update */
static double fluidstep = 0.0;

/***********************************************************************/

#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"

/********************** The Main LB Part *************************************/
void lb_init() {
  LB_TRACE(printf("Begin initializing fluid on CPU\n"));

  if (lbpar.agrid <= 0.0) {
    runtimeErrorMsg()
        << "Lattice Boltzmann agrid not set when initializing fluid";
  }
  if (check_runtime_errors(comm_cart))
    return;

  /* initialize the local lattice domain */

  try {
    lblattice = Lattice(lbpar.agrid, 0.5 /*offset*/, 1 /*halo size*/,
                        local_geo.length(), local_geo.my_right(),
                        box_geo.length(), calc_node_pos(comm_cart), node_grid);
  } catch (const std::runtime_error &e) {
    runtimeErrorMsg() << e.what();
    return;
  }

  /* allocate memory for data structures */
  lb_realloc_fluid();

  /* prepare the halo communication */
  lb_prepare_communication();

  /* initialize derived parameters */
  lb_reinit_parameters();

  /* setup the initial populations */
  lb_reinit_fluid();
}

void lb_reinit_fluid() {
  std::fill(lbfields.begin(), lbfields.end(), LB_FluidNode());
  /* default values for fields in lattice units */
  Utils::Vector3d momentum_density{};
  Utils::Vector6d stress{};

  for (Lattice::index_t index = 0; index < lblattice.halo_grid_volume;
       ++index) {
    // sets equilibrium distribution
    lb_set_population_from_density_momentum_density_stress(
        index, lbpar.density, momentum_density, stress);

#ifdef LB_BOUNDARIES
    lbfields[index].boundary = 0;
#endif // LB_BOUNDARIES
  }

#ifdef LB_BOUNDARIES
  LBBoundaries::lb_init_boundaries();
#endif // LB_BOUNDARIES
}

void lb_reinit_parameters() {
  if (lbpar.viscosity > 0.0) {
    /* Eq. (80) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar.gamma_shear = 1. - 2. / (6. * lbpar.viscosity + 1.);
  }

  if (lbpar.bulk_viscosity > 0.0) {
    /* Eq. (81) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar.gamma_bulk = 1. - 2. / (9. * lbpar.bulk_viscosity + 1.);
  }

  if (lbpar.is_TRT) {
    lbpar.gamma_bulk = lbpar.gamma_shear;
    lbpar.gamma_even = lbpar.gamma_shear;
    lbpar.gamma_odd =
        -(7.0 * lbpar.gamma_even + 1.0) / (lbpar.gamma_even + 7.0);
    // gamma_odd = lbpar.gamma_shear; //uncomment for BGK
  }

  // lbpar.gamma_shear = 0.0; //uncomment for special case of BGK
  // lbpar.gamma_bulk = 0.0;
  // gamma_odd = 0.0;
  // gamma_even = 0.0;

  if (lbpar.kT > 0.0) {
    /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).
     * Note that the modes are not normalized as in the paper here! */
    double mu = lbpar.kT / D3Q19::c_sound_sq<double> * lbpar.tau * lbpar.tau /
                (lbpar.agrid * lbpar.agrid);

    for (int i = 0; i < 4; i++)
      lbpar.phi[i] = 0.0;
    lbpar.phi[4] =
        sqrt(mu * D3Q19::w_k[4] * (1. - Utils::sqr(lbpar.gamma_bulk)));
    for (int i = 5; i < 10; i++)
      lbpar.phi[i] =
          sqrt(mu * D3Q19::w_k[i] * (1. - Utils::sqr(lbpar.gamma_shear)));
    for (int i = 10; i < 16; i++)
      lbpar.phi[i] =
          sqrt(mu * D3Q19::w_k[i] * (1 - Utils::sqr(lbpar.gamma_odd)));
    for (int i = 16; i < 19; i++)
      lbpar.phi[i] =
          sqrt(mu * D3Q19::w_k[i] * (1 - Utils::sqr(lbpar.gamma_even)));

    LB_TRACE(fprintf(
        stderr,
        "%d: lbpar.gamma_shear=%lf lbpar.gamma_bulk=%lf shear_fluct=%lf "
        "bulk_fluct=%lf mu=%lf, bulkvisc=%lf, visc=%lf\n",
        this_node, lbpar.gamma_shear, lbpar.gamma_bulk, lbpar.phi[9],
        lbpar.phi[4], mu, lbpar.bulk_viscosity, lbpar.viscosity));
  } else {
    for (int i = 0; i < D3Q19::n_vel; i++)
      lbpar.phi[i] = 0.0;
  }
}

/* Halo communication for push scheme */
static void halo_push_communication(LB_Fluid &lbfluid) {
  Lattice::index_t index;
  int x, y, z, count;
  int rnode, snode;
  double *buffer = nullptr, *sbuf = nullptr, *rbuf = nullptr;
  MPI_Status status;

  auto const yperiod = lblattice.halo_grid[0];
  auto const zperiod = lblattice.halo_grid[0] * lblattice.halo_grid[1];

  auto const node_neighbors = calc_node_neighbors(comm_cart);

  /***************
   * X direction *
   ***************/
  count = 5 * lblattice.halo_grid[1] * lblattice.halo_grid[2];
  sbuf = (double *)Utils::malloc(count * sizeof(double));
  rbuf = (double *)Utils::malloc(count * sizeof(double));

  /* send to right, recv from left i = 1, 7, 9, 11, 13 */
  snode = node_neighbors[1];
  rnode = node_neighbors[0];

  buffer = sbuf;
  index = get_linear_index(lblattice.grid[0] + 1, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      buffer[0] = lbfluid[1][index];
      buffer[1] = lbfluid[7][index];
      buffer[2] = lbfluid[9][index];
      buffer[3] = lbfluid[11][index];
      buffer[4] = lbfluid[13][index];
      buffer += 5;

      index += yperiod;
    }
  }

  MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
               MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);

  buffer = rbuf;
  index = get_linear_index(1, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      lbfluid[1][index] = buffer[0];
      lbfluid[7][index] = buffer[1];
      lbfluid[9][index] = buffer[2];
      lbfluid[11][index] = buffer[3];
      lbfluid[13][index] = buffer[4];
      buffer += 5;

      index += yperiod;
    }
  }

  /* send to left, recv from right i = 2, 8, 10, 12, 14 */
  snode = node_neighbors[0];
  rnode = node_neighbors[1];

  buffer = sbuf;
  index = get_linear_index(0, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      buffer[0] = lbfluid[2][index];
      buffer[1] = lbfluid[8][index];
      buffer[2] = lbfluid[10][index];
      buffer[3] = lbfluid[12][index];
      buffer[4] = lbfluid[14][index];
      buffer += 5;

      index += yperiod;
    }
  }

  MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
               MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);

  buffer = rbuf;
  index = get_linear_index(lblattice.grid[0], 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      lbfluid[2][index] = buffer[0];
      lbfluid[8][index] = buffer[1];
      lbfluid[10][index] = buffer[2];
      lbfluid[12][index] = buffer[3];
      lbfluid[14][index] = buffer[4];
      buffer += 5;

      index += yperiod;
    }
  }

  /***************
   * Y direction *
   ***************/
  count = 5 * lblattice.halo_grid[0] * lblattice.halo_grid[2];
  sbuf = Utils::realloc(sbuf, count * sizeof(double));
  rbuf = Utils::realloc(rbuf, count * sizeof(double));

  /* send to right, recv from left i = 3, 7, 10, 15, 17 */
  snode = node_neighbors[3];
  rnode = node_neighbors[2];

  buffer = sbuf;
  index = get_linear_index(0, lblattice.grid[1] + 1, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[3][index];
      buffer[1] = lbfluid[7][index];
      buffer[2] = lbfluid[10][index];
      buffer[3] = lbfluid[15][index];
      buffer[4] = lbfluid[17][index];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
               MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);

  buffer = rbuf;
  index = get_linear_index(0, 1, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[3][index] = buffer[0];
      lbfluid[7][index] = buffer[1];
      lbfluid[10][index] = buffer[2];
      lbfluid[15][index] = buffer[3];
      lbfluid[17][index] = buffer[4];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  /* send to left, recv from right i = 4, 8, 9, 16, 18 */
  snode = node_neighbors[2];
  rnode = node_neighbors[3];

  buffer = sbuf;
  index = get_linear_index(0, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[4][index];
      buffer[1] = lbfluid[8][index];
      buffer[2] = lbfluid[9][index];
      buffer[3] = lbfluid[16][index];
      buffer[4] = lbfluid[18][index];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
               MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);

  buffer = rbuf;
  index = get_linear_index(0, lblattice.grid[1], 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[4][index] = buffer[0];
      lbfluid[8][index] = buffer[1];
      lbfluid[9][index] = buffer[2];
      lbfluid[16][index] = buffer[3];
      lbfluid[18][index] = buffer[4];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  /***************
   * Z direction *
   ***************/
  count = 5 * lblattice.halo_grid[0] * lblattice.halo_grid[1];
  sbuf = Utils::realloc(sbuf, count * sizeof(double));
  rbuf = Utils::realloc(rbuf, count * sizeof(double));

  /* send to right, recv from left i = 5, 11, 14, 15, 18 */
  snode = node_neighbors[5];
  rnode = node_neighbors[4];

  buffer = sbuf;
  index = get_linear_index(0, 0, lblattice.grid[2] + 1, lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[5][index];
      buffer[1] = lbfluid[11][index];
      buffer[2] = lbfluid[14][index];
      buffer[3] = lbfluid[15][index];
      buffer[4] = lbfluid[18][index];
      buffer += 5;

      ++index;
    }
  }

  MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
               MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);

  buffer = rbuf;
  index = get_linear_index(0, 0, 1, lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[5][index] = buffer[0];
      lbfluid[11][index] = buffer[1];
      lbfluid[14][index] = buffer[2];
      lbfluid[15][index] = buffer[3];
      lbfluid[18][index] = buffer[4];
      buffer += 5;

      ++index;
    }
  }

  /* send to left, recv from right i = 6, 12, 13, 16, 17 */
  snode = node_neighbors[4];
  rnode = node_neighbors[5];

  buffer = sbuf;
  index = get_linear_index(0, 0, 0, lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[6][index];
      buffer[1] = lbfluid[12][index];
      buffer[2] = lbfluid[13][index];
      buffer[3] = lbfluid[16][index];
      buffer[4] = lbfluid[17][index];
      buffer += 5;

      ++index;
    }
  }

  MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
               MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);

  buffer = rbuf;
  index = get_linear_index(0, 0, lblattice.grid[2], lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[6][index] = buffer[0];
      lbfluid[12][index] = buffer[1];
      lbfluid[13][index] = buffer[2];
      lbfluid[16][index] = buffer[3];
      lbfluid[17][index] = buffer[4];
      buffer += 5;

      ++index;
    }
  }

  free(rbuf);
  free(sbuf);
}

/***********************************************************************/

/** Performs basic sanity checks. */
void lb_sanity_checks() {

  if (lbpar.agrid <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann agrid not set";
  }
  if (lbpar.tau <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann time step not set";
  }
  if (lbpar.density <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann fluid density not set";
  }
  if (lbpar.viscosity <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann fluid viscosity not set";
  }
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    runtimeErrorMsg() << "LB requires domain-decomposition cellsystem";
  }
  if (skin == 0.0) {
    runtimeErrorMsg() << "LB requires a positive skin";
  }
  if (cell_structure.use_verlet_list && skin >= lbpar.agrid / 2.0) {
    runtimeErrorMsg() << "LB requires either no Verlet lists or that the skin "
                         "of the verlet list to be less than half of "
                         "lattice-Boltzmann grid spacing";
  }
}

uint64_t lb_fluid_get_rng_state() {
  assert(rng_counter_fluid);
  return rng_counter_fluid->value();
}

void mpi_set_lb_fluid_counter(uint64_t counter) {
  rng_counter_fluid = Utils::Counter<uint64_t>(counter);
}

REGISTER_CALLBACK(mpi_set_lb_fluid_counter)

void lb_fluid_set_rng_state(uint64_t counter) {
  mpi_call(mpi_set_lb_fluid_counter, counter);
  mpi_set_lb_fluid_counter(counter);
}

/***********************************************************************/

/** (Re-)allocate memory for the fluid and initialize pointers. */
void lb_realloc_fluid() {
  LB_TRACE(printf("reallocating fluid\n"));
  const std::array<int, 2> size = {{D3Q19::n_vel, lblattice.halo_grid_volume}};

  lbfluid_a.resize(size);
  lbfluid_b.resize(size);

  using Utils::Span;
  for (int i = 0; i < size[0]; i++) {
    lbfluid[i] = Span<double>(lbfluid_a[i].origin(), size[1]);
    lbfluid_post[i] = Span<double>(lbfluid_b[i].origin(), size[1]);
  }

  lbfields.resize(lblattice.halo_grid_volume);
}

/** Set up the structures for exchange of the halo regions.
 *  See also \ref halo.cpp
 */
void lb_prepare_communication() {
  int i;
  HaloCommunicator comm = {0, nullptr};

  /* since the data layout is a structure of arrays, we have to
   * generate a communication for this structure: first we generate
   * the communication for one of the arrays (the 0-th velocity
   * population), then we replicate this communication for the other
   * velocity indices by constructing appropriate vector
   * datatypes */

  /* prepare the communication for a single velocity */
  prepare_halo_communication(&comm, &lblattice, FIELDTYPE_DOUBLE, MPI_DOUBLE,
                             node_grid);

  update_halo_comm.num = comm.num;
  update_halo_comm.halo_info =
      Utils::realloc(update_halo_comm.halo_info, comm.num * sizeof(HaloInfo));

  /* replicate the halo structure */
  for (i = 0; i < comm.num; i++) {
    HaloInfo *hinfo = &(update_halo_comm.halo_info[i]);

    hinfo->source_node = comm.halo_info[i].source_node;
    hinfo->dest_node = comm.halo_info[i].dest_node;
    hinfo->s_offset = comm.halo_info[i].s_offset;
    hinfo->r_offset = comm.halo_info[i].r_offset;
    hinfo->type = comm.halo_info[i].type;

    /* generate the vector datatype for the structure of lattices we
     * have to use hvector here because the extent of the subtypes
     * does not span the full lattice and hence we cannot get the
     * correct vskip out of them */

    MPI_Aint lower;
    MPI_Aint extent;
    MPI_Type_get_extent(MPI_DOUBLE, &lower, &extent);
    MPI_Type_create_hvector(D3Q19::n_vel, 1,
                            lblattice.halo_grid_volume * extent,
                            comm.halo_info[i].datatype, &hinfo->datatype);
    MPI_Type_commit(&hinfo->datatype);

    halo_create_field_hvector(D3Q19::n_vel, 1,
                              lblattice.halo_grid_volume * sizeof(double),
                              comm.halo_info[i].fieldtype, &hinfo->fieldtype);
  }

  release_halo_communication(&comm);
}

/***********************************************************************/
/** \name Mapping between hydrodynamic fields and particle populations */
/***********************************************************************/
/*@{*/
Utils::Vector19d lb_get_population_from_density_momentum_density_stress(
    double density, Utils::Vector3d const &momentum_density,
    Utils::Vector6d const &stress) {
  Utils::Vector19d population{};
  auto const trace = stress[0] + stress[2] + stress[5];

  /* update the q=0 sublattice */
  population[0] = 1. / 3. * density - 1. / 2. * trace;

  /* update the q=1 sublattice */
  auto density_times_coeff = 1. / 18. * density;

  population[1] = density_times_coeff + 1. / 6. * momentum_density[0] +
                  1. / 4. * stress[0] - 1. / 12. * trace;
  population[2] = density_times_coeff - 1. / 6. * momentum_density[0] +
                  1. / 4. * stress[0] - 1. / 12. * trace;
  population[3] = density_times_coeff + 1. / 6. * momentum_density[1] +
                  1. / 4. * stress[2] - 1. / 12. * trace;
  population[4] = density_times_coeff - 1. / 6. * momentum_density[1] +
                  1. / 4. * stress[2] - 1. / 12. * trace;
  population[5] = density_times_coeff + 1. / 6. * momentum_density[2] +
                  1. / 4. * stress[5] - 1. / 12. * trace;
  population[6] = density_times_coeff - 1. / 6. * momentum_density[2] +
                  1. / 4. * stress[5] - 1. / 12. * trace;

  /* update the q=2 sublattice */
  density_times_coeff = 1. / 36. * density;

  auto tmp1 = stress[0] + stress[2];
  auto tmp2 = 2.0 * stress[1];

  population[7] = density_times_coeff +
                  1. / 12. * (momentum_density[0] + momentum_density[1]) +
                  1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  population[8] = density_times_coeff -
                  1. / 12. * (momentum_density[0] + momentum_density[1]) +
                  1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  population[9] = density_times_coeff +
                  1. / 12. * (momentum_density[0] - momentum_density[1]) +
                  1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  population[10] = density_times_coeff -
                   1. / 12. * (momentum_density[0] - momentum_density[1]) +
                   1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;

  tmp1 = stress[0] + stress[5];
  tmp2 = 2.0 * stress[3];

  population[11] = density_times_coeff +
                   1. / 12. * (momentum_density[0] + momentum_density[2]) +
                   1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  population[12] = density_times_coeff -
                   1. / 12. * (momentum_density[0] + momentum_density[2]) +
                   1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  population[13] = density_times_coeff +
                   1. / 12. * (momentum_density[0] - momentum_density[2]) +
                   1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  population[14] = density_times_coeff -
                   1. / 12. * (momentum_density[0] - momentum_density[2]) +
                   1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;

  tmp1 = stress[2] + stress[5];
  tmp2 = 2.0 * stress[4];

  population[15] = density_times_coeff +
                   1. / 12. * (momentum_density[1] + momentum_density[2]) +
                   1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  population[16] = density_times_coeff -
                   1. / 12. * (momentum_density[1] + momentum_density[2]) +
                   1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  population[17] = density_times_coeff +
                   1. / 12. * (momentum_density[1] - momentum_density[2]) +
                   1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  population[18] = density_times_coeff -
                   1. / 12. * (momentum_density[1] - momentum_density[2]) +
                   1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  return population;
}

void lb_set_population_from_density_momentum_density_stress(
    Lattice::index_t const index, double density,
    Utils::Vector3d const &momentum_density, Utils::Vector6d const &stress) {
  auto const population =
      lb_get_population_from_density_momentum_density_stress(
          density, momentum_density, stress);
  lb_set_population(index, population);
}
/*@}*/

#include <boost/range/numeric.hpp>

/** Calculation of hydrodynamic modes */
std::array<double, 19> lb_calc_modes(Lattice::index_t index) {
  return Utils::matrix_vector_product<double, 19, e_ki>(
      LB_Fluid_Ref(index, lbfluid));
}

template <typename T>
inline std::array<T, 19> lb_relax_modes(Lattice::index_t index,
                                        const std::array<T, 19> &modes) {
  T density, momentum_density[3], stress_eq[6];

  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */
  density = modes[0] + lbpar.density;

  momentum_density[0] = modes[1] + 0.5 * lbfields[index].force_density[0];
  momentum_density[1] = modes[2] + 0.5 * lbfields[index].force_density[1];
  momentum_density[2] = modes[3] + 0.5 * lbfields[index].force_density[2];

  using Utils::sqr;
  auto const momentum_density2 = sqr(momentum_density[0]) +
                                 sqr(momentum_density[1]) +
                                 sqr(momentum_density[2]);

  /* equilibrium part of the stress modes */
  stress_eq[0] = momentum_density2 / density;
  stress_eq[1] =
      (sqr(momentum_density[0]) - sqr(momentum_density[1])) / density;
  stress_eq[2] = (momentum_density2 - 3.0 * sqr(momentum_density[2])) / density;
  stress_eq[3] = momentum_density[0] * momentum_density[1] / density;
  stress_eq[4] = momentum_density[0] * momentum_density[2] / density;
  stress_eq[5] = momentum_density[1] * momentum_density[2] / density;

  return {{modes[0], modes[1], modes[2], modes[3],
           /* relax the stress modes */
           stress_eq[0] + lbpar.gamma_bulk * (modes[4] - stress_eq[0]),
           stress_eq[1] + lbpar.gamma_shear * (modes[5] - stress_eq[1]),
           stress_eq[2] + lbpar.gamma_shear * (modes[6] - stress_eq[2]),
           stress_eq[3] + lbpar.gamma_shear * (modes[7] - stress_eq[3]),
           stress_eq[4] + lbpar.gamma_shear * (modes[8] - stress_eq[4]),
           stress_eq[5] + lbpar.gamma_shear * (modes[9] - stress_eq[5]),
           /* relax the ghost modes (project them out) */
           /* ghost modes have no equilibrium part due to orthogonality */
           lbpar.gamma_odd * modes[10], lbpar.gamma_odd * modes[11],
           lbpar.gamma_odd * modes[12], lbpar.gamma_odd * modes[13],
           lbpar.gamma_odd * modes[14], lbpar.gamma_odd * modes[15],
           lbpar.gamma_even * modes[16], lbpar.gamma_even * modes[17],
           lbpar.gamma_even * modes[18]}};
}

template <typename T>
inline std::array<T, 19> lb_thermalize_modes(Lattice::index_t index,
                                             const std::array<T, 19> &modes) {
  if (lbpar.kT > 0.0) {
    using Utils::uniform;
    using rng_type = r123::Philox4x64;
    using ctr_type = rng_type::ctr_type;

    const r123::Philox4x64::ctr_type c{
        {rng_counter_fluid->value(), static_cast<uint64_t>(RNGSalt::FLUID)}};
    const T rootdensity = std::sqrt(std::fabs(modes[0] + lbpar.density));
    auto const pref = std::sqrt(12.) * rootdensity;

    const ctr_type noise[4] = {
        rng_type{}(c, {{static_cast<uint64_t>(index), 0ul}}),
        rng_type{}(c, {{static_cast<uint64_t>(index), 1ul}}),
        rng_type{}(c, {{static_cast<uint64_t>(index), 2ul}}),
        rng_type{}(c, {{static_cast<uint64_t>(index), 3ul}})};

    auto rng = [&](int i) { return uniform(noise[i / 4][i % 4]); };

    return {/* conserved modes */
            {modes[0], modes[1], modes[2], modes[3],
             /* stress modes */
             modes[4] + pref * lbpar.phi[4] * rng(0),
             modes[5] + pref * lbpar.phi[5] * rng(1),
             modes[6] + pref * lbpar.phi[6] * rng(2),
             modes[7] + pref * lbpar.phi[7] * rng(3),
             modes[8] + pref * lbpar.phi[8] * rng(4),
             modes[9] + pref * lbpar.phi[9] * rng(5),

             /* ghost modes */
             modes[10] + pref * lbpar.phi[10] * rng(6),
             modes[11] + pref * lbpar.phi[11] * rng(7),
             modes[12] + pref * lbpar.phi[12] * rng(8),
             modes[13] + pref * lbpar.phi[13] * rng(9),
             modes[14] + pref * lbpar.phi[14] * rng(10),
             modes[15] + pref * lbpar.phi[15] * rng(11),
             modes[16] + pref * lbpar.phi[16] * rng(12),
             modes[17] + pref * lbpar.phi[17] * rng(13),
             modes[18] + pref * lbpar.phi[18] * rng(14)}};
  }
  return modes;
}

template <typename T>
std::array<T, 19> lb_apply_forces(Lattice::index_t index,
                                  const std::array<T, 19> &modes) {
  const auto &f = lbfields[index].force_density;

  auto const density = modes[0] + lbpar.density;

  /* hydrodynamic momentum density is redefined when external forces present */
  auto const u = Utils::Vector3d{modes[1] + 0.5 * f[0], modes[2] + 0.5 * f[1],
                                 modes[3] + 0.5 * f[2]} /
                 density;

  double C[6];
  C[0] = (1. + lbpar.gamma_bulk) * u[0] * f[0] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * (u * f);
  C[2] = (1. + lbpar.gamma_bulk) * u[1] * f[1] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * (u * f);
  C[5] = (1. + lbpar.gamma_bulk) * u[2] * f[2] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * (u * f);
  C[1] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[0] * f[1] + u[1] * f[0]);
  C[3] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[0] * f[2] + u[2] * f[0]);
  C[4] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[1] * f[2] + u[2] * f[1]);

  return {{modes[0],
           /* update momentum modes */
           modes[1] + f[0], modes[2] + f[1], modes[3] + f[2],
           /* update stress modes */
           modes[4] + C[0] + C[2] + C[5], modes[5] + C[0] - C[2],
           modes[6] + C[0] + C[2] - 2. * C[5], modes[7] + C[1], modes[8] + C[3],
           modes[9] + C[4], modes[10], modes[11], modes[12], modes[13],
           modes[14], modes[15], modes[16], modes[17], modes[18]}};
}

template <typename T>
std::array<T, 19> normalize_modes(const std::array<T, 19> &modes) {
  auto normalized_modes = modes;
  for (int i = 0; i < modes.size(); i++) {
    normalized_modes[i] /= D3Q19::w_k[i];
  }
  return normalized_modes;
}

template <typename T, std::size_t N>
std::array<T, N> lb_calc_n_from_m(const std::array<T, N> &modes) {
  auto ret = Utils::matrix_vector_product<T, N, e_ki_transposed>(
      normalize_modes(modes));
  std::transform(ret.begin(), ret.end(), ::D3Q19::w.begin(), ret.begin(),
                 std::multiplies<T>());
  return ret;
}

template <typename T>
inline void lb_calc_n_from_modes_push(LB_Fluid &lbfluid, Lattice::index_t index,
                                      const std::array<T, 19> &m) {
  const std::array<int, 3> period = {
      {1, lblattice.halo_grid[0],
       lblattice.halo_grid[0] * lblattice.halo_grid[1]}};
  auto const f = lb_calc_n_from_m(m);
  for (int i = 0; i < 19; i++) {
    auto const next = index + boost::inner_product(period, D3Q19::c[i], 0);
    lbfluid[i][next] = f[i];
  }
}

/* Collisions and streaming (push scheme) */
inline void lb_collide_stream() {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
/* loop over all lattice cells (halo excluded) */
#ifdef LB_BOUNDARIES
  for (auto &lbboundarie : LBBoundaries::lbboundaries) {
    (*lbboundarie).reset_force();
  }
#endif // LB_BOUNDARIES

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  // Safeguard the node forces so that we can later use them for the IBM
  // particle update
  // In the following loop the lbfields[XX].force are reset to zero
  // Safeguard the node forces so that we can later use them for the IBM
  // particle update In the following loop the lbfields[XX].force are reset to
  // zero
  for (int i = 0; i < lblattice.halo_grid_volume; ++i) {
    lbfields[i].force_density_buf = lbfields[i].force_density;
  }
#endif

  Lattice::index_t index = lblattice.halo_offset;
  for (int z = 1; z <= lblattice.grid[2]; z++) {
    for (int y = 1; y <= lblattice.grid[1]; y++) {
      for (int x = 1; x <= lblattice.grid[0]; x++) {
// as we only want to apply this to non-boundary nodes we can throw out
// the if-clause if we have a non-bounded domain
#ifdef LB_BOUNDARIES
        if (!lbfields[index].boundary)
#endif // LB_BOUNDARIES
        {
          /* calculate modes locally */
          auto const modes = lb_calc_modes(index);

          /* deterministic collisions */
          auto const relaxed_modes = lb_relax_modes(index, modes);

          /* fluctuating hydrodynamics */
          auto const thermalized_modes =
              lb_thermalize_modes(index, relaxed_modes);

          /* apply forces */
          auto const modes_with_forces =
              lb_apply_forces(index, thermalized_modes);

          /* reset the force density */
          lbfields[index].force_density = lbpar.ext_force_density;

          /* transform back to populations and streaming */
          lb_calc_n_from_modes_push(lbfluid_post, index, modes_with_forces);
        }

        ++index; /* next node */
      }
      index += 2; /* skip halo region */
    }
    index += 2 * lblattice.halo_grid[0]; /* skip halo region */
  }

  /* exchange halo regions */
  halo_push_communication(lbfluid_post);

#ifdef LB_BOUNDARIES
  /* boundary conditions for links */
  lb_bounce_back(lbfluid_post);
#endif // LB_BOUNDARIES

  /* swap the pointers for old and new population fields */
  std::swap(lbfluid, lbfluid_post);

  halo_communication(&update_halo_comm,
                     reinterpret_cast<char *>(lbfluid[0].data()));

#ifdef ADDITIONAL_CHECKS
  lb_check_halo_regions(lbfluid);
#endif
}

/***********************************************************************/
/** \name Update step for the lattice Boltzmann fluid                  */
/***********************************************************************/
/*@{*/

/** Update the lattice Boltzmann fluid.
 *
 *  This function is called from the integrator. Since the time step
 *  for the lattice dynamics can be coarser than the MD time step, we
 *  monitor the time since the last lattice update.
 */
void lattice_boltzmann_update() {
  auto factor = (int)round(lbpar.tau / time_step);

  fluidstep += 1;
  if (fluidstep >= factor) {
    fluidstep = 0;

    lb_collide_stream();
  }
}

/*@}*/

/***********************************************************************/
/** \name Coupling part */
/***********************************************************************/
/*@{*/

/***********************************************************************/

static int compare_buffers(double *buf1, double *buf2, int size) {
  int ret;
  if (memcmp(buf1, buf2, size) != 0) {
    runtimeErrorMsg() << "Halo buffers are not identical";
    ret = 1;
  } else {
    ret = 0;
  }
  return ret;
}

/** Check consistency of the halo regions (ADDITIONAL_CHECKS)
 *  This function can be used as an additional check. It test whether the
 *  halo regions have been exchanged correctly.
 */
void lb_check_halo_regions(const LB_Fluid &lbfluid) {
  Lattice::index_t index;
  int i, x, y, z, s_node, r_node, count = D3Q19::n_vel;
  double *s_buffer, *r_buffer;
  MPI_Status status[2];

  r_buffer = (double *)Utils::malloc(count * sizeof(double));
  s_buffer = (double *)Utils::malloc(count * sizeof(double));

  auto const node_neighbors = calc_node_neighbors(comm_cart);

  if (box_geo.periodic(0)) {
    for (z = 0; z < lblattice.halo_grid[2]; ++z) {
      for (y = 0; y < lblattice.halo_grid[1]; ++y) {
        index = get_linear_index(0, y, z, lblattice.halo_grid);
        for (i = 0; i < D3Q19::n_vel; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[1];
        r_node = node_neighbors[0];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index =
              get_linear_index(lblattice.grid[0], y, z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index =
              get_linear_index(lblattice.grid[0], y, z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=" << 0 << " at index=" << index
                      << " y=" << y << " z=" << z << "\n";
          }
        }

        index =
            get_linear_index(lblattice.grid[0] + 1, y, z, lblattice.halo_grid);
        for (i = 0; i < D3Q19::n_vel; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[0];
        r_node = node_neighbors[1];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index = get_linear_index(1, y, z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index = get_linear_index(1, y, z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=0 at index=" << index
                      << " y=" << y << " z=" << z << "\n";
          }
        }
      }
    }
  }

  if (box_geo.periodic(1)) {
    for (z = 0; z < lblattice.halo_grid[2]; ++z) {
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index = get_linear_index(x, 0, z, lblattice.halo_grid);
        for (i = 0; i < D3Q19::n_vel; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[3];
        r_node = node_neighbors[2];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index =
              get_linear_index(x, lblattice.grid[1], z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index =
              get_linear_index(x, lblattice.grid[1], z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=1 at index=" << index
                      << " x=" << x << " z=" << z << "\n";
          }
        }
      }
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index =
            get_linear_index(x, lblattice.grid[1] + 1, z, lblattice.halo_grid);
        for (i = 0; i < D3Q19::n_vel; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[2];
        r_node = node_neighbors[3];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index = get_linear_index(x, 1, z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index = get_linear_index(x, 1, z, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=1 at index=" << index
                      << " x=" << x << " z=" << z << "\n";
          }
        }
      }
    }
  }

  if (box_geo.periodic(2)) {
    for (y = 0; y < lblattice.halo_grid[1]; ++y) {
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index = get_linear_index(x, y, 0, lblattice.halo_grid);
        for (i = 0; i < D3Q19::n_vel; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[5];
        r_node = node_neighbors[4];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index =
              get_linear_index(x, y, lblattice.grid[2], lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index =
              get_linear_index(x, y, lblattice.grid[2], lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=2 at index=" << index
                      << " x=" << x << " y=" << y << " z=" << lblattice.grid[2]
                      << "\n";
          }
        }
      }
    }
    for (y = 0; y < lblattice.halo_grid[1]; ++y) {
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index =
            get_linear_index(x, y, lblattice.grid[2] + 1, lblattice.halo_grid);
        for (i = 0; i < D3Q19::n_vel; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[4];
        r_node = node_neighbors[5];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index = get_linear_index(x, y, 1, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index = get_linear_index(x, y, 1, lblattice.halo_grid);
          for (i = 0; i < D3Q19::n_vel; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=2 at index=" << index
                      << " x=" << x << " y=" << y << "\n";
          }
        }
      }
    }
  }

  free(r_buffer);
  free(s_buffer);
}

double lb_calc_density(std::array<double, 19> const &modes) {
  return modes[0] + lbpar.density;
}

Utils::Vector3d lb_calc_momentum_density(std::array<double, 19> const &modes,
                                         Utils::Vector3d const &force_density) {
  return Utils::Vector3d{{modes[1] + 0.5 * force_density[0],
                          modes[2] + 0.5 * force_density[1],
                          modes[3] + 0.5 * force_density[2]}};
}

Utils::Vector6d lb_calc_stress(std::array<double, 19> const &modes,
                               Utils::Vector3d const &force_density) {
  auto const momentum_density = lb_calc_momentum_density(modes, force_density);
  auto const density = lb_calc_density(modes);
  using Utils::sqr;
  auto const momentum_density2 = sqr(momentum_density[0]) +
                                 sqr(momentum_density[1]) +
                                 sqr(momentum_density[2]);
  /* equilibrium part of the stress modes */
  Utils::Vector6d modes_from_stress_eq{};
  modes_from_stress_eq[0] = momentum_density2 / density;
  modes_from_stress_eq[1] =
      (sqr(momentum_density[0]) - sqr(momentum_density[1])) / density;
  modes_from_stress_eq[2] =
      (momentum_density2 - 3.0 * sqr(momentum_density[2])) / density;
  modes_from_stress_eq[3] = momentum_density[0] * momentum_density[1] / density;
  modes_from_stress_eq[4] = momentum_density[0] * momentum_density[2] / density;
  modes_from_stress_eq[5] = momentum_density[1] * momentum_density[2] / density;

  /* Now we must predict the outcome of the next collision */
  /* We immediately average pre- and post-collision. */

  Utils::Vector6d avg_modes;
  avg_modes[0] =
      modes_from_stress_eq[0] +
      (0.5 + 0.5 * lbpar.gamma_bulk) * (modes[4] - modes_from_stress_eq[0]);
  avg_modes[1] =
      modes_from_stress_eq[1] +
      (0.5 + 0.5 * lbpar.gamma_shear) * (modes[5] - modes_from_stress_eq[1]);
  avg_modes[2] =
      modes_from_stress_eq[2] +
      (0.5 + 0.5 * lbpar.gamma_shear) * (modes[6] - modes_from_stress_eq[2]);
  avg_modes[3] =
      modes_from_stress_eq[3] +
      (0.5 + 0.5 * lbpar.gamma_shear) * (modes[7] - modes_from_stress_eq[3]);
  avg_modes[4] =
      modes_from_stress_eq[4] +
      (0.5 + 0.5 * lbpar.gamma_shear) * (modes[8] - modes_from_stress_eq[4]);
  avg_modes[5] =
      modes_from_stress_eq[5] +
      (0.5 + 0.5 * lbpar.gamma_shear) * (modes[9] - modes_from_stress_eq[5]);

  // Transform the stress tensor components according to the modes that
  // correspond to those used by U. Schiller. In terms of populations this
  // expression then corresponds exactly to those in Eqs. 116 - 121 in the
  // Duenweg and Ladd paper, when these are written out in populations.
  // But to ensure this, the expression in Schiller's modes has to be different!

  Utils::Vector6d stress;
  stress[0] =
      (2.0 * (modes[0] + avg_modes[0]) + avg_modes[2] + 3.0 * avg_modes[1]) /
      6.0;                  // xx
  stress[1] = avg_modes[3]; // xy
  stress[2] =
      (2.0 * (modes[0] + avg_modes[0]) + avg_modes[2] - 3.0 * avg_modes[1]) /
      6.0;                                                    // yy
  stress[3] = avg_modes[4];                                   // xz
  stress[4] = avg_modes[5];                                   // yz
  stress[5] = (modes[0] + avg_modes[0] - avg_modes[2]) / 3.0; // zz
  return stress;
}

#ifdef LB_BOUNDARIES
/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 */
void lb_bounce_back(LB_Fluid &lbfluid) {
  int k, i, l;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0] * lblattice.halo_grid[1];
  int next[19];
  double population_shift;
  next[0] = 0;                     // ( 0, 0, 0) =
  next[1] = 1;                     // ( 1, 0, 0) +
  next[2] = -1;                    // (-1, 0, 0)
  next[3] = yperiod;               // ( 0, 1, 0) +
  next[4] = -yperiod;              // ( 0,-1, 0)
  next[5] = zperiod;               // ( 0, 0, 1) +
  next[6] = -zperiod;              // ( 0, 0,-1)
  next[7] = (1 + yperiod);         // ( 1, 1, 0) +
  next[8] = -(1 + yperiod);        // (-1,-1, 0)
  next[9] = (1 - yperiod);         // ( 1,-1, 0)
  next[10] = -(1 - yperiod);       // (-1, 1, 0) +
  next[11] = (1 + zperiod);        // ( 1, 0, 1) +
  next[12] = -(1 + zperiod);       // (-1, 0,-1)
  next[13] = (1 - zperiod);        // ( 1, 0,-1)
  next[14] = -(1 - zperiod);       // (-1, 0, 1) +
  next[15] = (yperiod + zperiod);  // ( 0, 1, 1) +
  next[16] = -(yperiod + zperiod); // ( 0,-1,-1)
  next[17] = (yperiod - zperiod);  // ( 0, 1,-1)
  next[18] = -(yperiod - zperiod); // ( 0,-1, 1) +
  int reverse[] = {0, 2,  1,  4,  3,  6,  5,  8,  7, 10,
                   9, 12, 11, 14, 13, 16, 15, 18, 17};

  /* bottom-up sweep */
  for (int z = 0; z < lblattice.grid[2] + 2; z++) {
    for (int y = 0; y < lblattice.grid[1] + 2; y++) {
      for (int x = 0; x < lblattice.grid[0] + 2; x++) {
        k = get_linear_index(x, y, z, lblattice.halo_grid);

        if (lbfields[k].boundary) {
          for (i = 0; i < 19; i++) {
            population_shift = 0;
            for (l = 0; l < 3; l++) {
              population_shift -= lbpar.density * 2 * D3Q19::c[i][l] *
                                  D3Q19::w[i] * lbfields[k].slip_velocity[l] /
                                  D3Q19::c_sound_sq<double>;
            }

            if (x - D3Q19::c[i][0] > 0 &&
                x - D3Q19::c[i][0] < lblattice.grid[0] + 1 &&
                y - D3Q19::c[i][1] > 0 &&
                y - D3Q19::c[i][1] < lblattice.grid[1] + 1 &&
                z - D3Q19::c[i][2] > 0 &&
                z - D3Q19::c[i][2] < lblattice.grid[2] + 1) {
              if (!lbfields[k - next[i]].boundary) {
                for (l = 0; l < 3; l++) {
                  (*LBBoundaries::lbboundaries[lbfields[k].boundary - 1])
                      .force()[l] += // TODO
                      (2 * lbfluid[i][k] + population_shift) * D3Q19::c[i][l];
                }
                lbfluid[reverse[i]][k - next[i]] =
                    lbfluid[i][k] + population_shift;
              } else {
                lbfluid[reverse[i]][k - next[i]] = lbfluid[i][k] = 0.0;
              }
            }
          }
        }
      }
    }
  }
}
#endif

/** Calculate the local fluid momentum.
 *  The calculation is implemented explicitly for the special case of D3Q19.
 *  @param[in]  index  Local lattice site
 *  @retval The local fluid momentum.
 */
inline Utils::Vector3d lb_calc_local_momentum_density(Lattice::index_t index) {
  return {{lbfluid[1][index] - lbfluid[2][index] + lbfluid[7][index] -
               lbfluid[8][index] + lbfluid[9][index] - lbfluid[10][index] +
               lbfluid[11][index] - lbfluid[12][index] + lbfluid[13][index] -
               lbfluid[14][index],
           lbfluid[3][index] - lbfluid[4][index] + lbfluid[7][index] -
               lbfluid[8][index] - lbfluid[9][index] + lbfluid[10][index] +
               lbfluid[15][index] - lbfluid[16][index] + lbfluid[17][index] -
               lbfluid[18][index],
           lbfluid[5][index] - lbfluid[6][index] + lbfluid[11][index] -
               lbfluid[12][index] - lbfluid[13][index] + lbfluid[14][index] +
               lbfluid[15][index] - lbfluid[16][index] - lbfluid[17][index] +
               lbfluid[18][index]}};
}

// Statistics in MD units.

/** Calculate momentum of the LB fluid.
 * \param result Fluid momentum
 */
void lb_calc_fluid_momentum(double *result) {

  int x, y, z, index;
  Utils::Vector3d momentum_density{}, momentum{};

  for (x = 1; x <= lblattice.grid[0]; x++) {
    for (y = 1; y <= lblattice.grid[1]; y++) {
      for (z = 1; z <= lblattice.grid[2]; z++) {
        index = get_linear_index(x, y, z, lblattice.halo_grid);

        momentum_density = lb_calc_local_momentum_density(index);
        momentum += momentum_density + .5 * lbfields[index].force_density;
      }
    }
  }

  momentum *= lbpar.agrid / lbpar.tau;

  MPI_Reduce(momentum.data(), result, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

void lb_collect_boundary_forces(double *result) {
#ifdef LB_BOUNDARIES
  int n_lb_boundaries = LBBoundaries::lbboundaries.size();
  std::vector<double> boundary_forces(3 * n_lb_boundaries);
  int i = 0;
  for (auto it = LBBoundaries::lbboundaries.begin();
       it != LBBoundaries::lbboundaries.end(); ++it, i++)
    for (int j = 0; j < 3; j++)
      boundary_forces[3 * i + j] = (**it).force()[j];

  MPI_Reduce(boundary_forces.data(), result, 3 * n_lb_boundaries, MPI_DOUBLE,
             MPI_SUM, 0, comm_cart);
#endif
}

/*@}*/
