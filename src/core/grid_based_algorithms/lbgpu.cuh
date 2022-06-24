/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
 *  Implementation in lbgpu_cuda.cu.
 */

#ifndef LBGPU_CUH
#define LBGPU_CUH

#include "config.hpp"

#ifdef CUDA
#include <curand_kernel.h>

#include <utils/Array.hpp>

/** Velocity densities for the lattice Boltzmann system. */
struct LB_nodes_gpu {
  /** velocity density of the node */
  Utils::Array<float, 19> *populations = nullptr;
  unsigned int *boundary = nullptr;
  Utils::Array<float, 3> *boundary_velocity = nullptr;
};

__device__ __inline__ float
calc_mode_x_from_n(Utils::Array<float, 19> const &populations, int x) {
  switch (x) {
  case 0:
    return populations[0] + populations[1] + populations[2] + populations[3] +
           populations[4] + populations[5] + populations[6] + populations[7] +
           populations[8] + populations[9] + populations[10] + populations[11] +
           populations[12] + populations[13] + populations[14] +
           populations[15] + populations[16] + populations[17] +
           populations[18];
  case 1:
    return (populations[1] - populations[2]) +
           (populations[7] - populations[8]) +
           (populations[9] - populations[10]) +
           (populations[11] - populations[12]) +
           (populations[13] - populations[14]);
  case 2:
    return (populations[3] - populations[4]) +
           (populations[7] - populations[8]) -
           (populations[9] - populations[10]) +
           (populations[15] - populations[16]) +
           (populations[17] - populations[18]);
  case 3:
    return (populations[5] - populations[6]) +
           (populations[11] - populations[12]) -
           (populations[13] - populations[14]) +
           (populations[15] - populations[16]) -
           (populations[17] - populations[18]);
  case 4:
    return -populations[0] + populations[7] + populations[8] + populations[9] +
           populations[10] + populations[11] + populations[12] +
           populations[13] + populations[14] + populations[15] +
           populations[16] + populations[17] + populations[18];
  case 5:
    return (populations[1] + populations[2]) -
           (populations[3] + populations[4]) +
           (populations[11] + populations[12]) +
           (populations[13] + populations[14]) -
           (populations[15] + populations[16]) -
           (populations[17] + populations[18]);
  case 6:
    return (populations[1] + populations[2]) +
           (populations[3] + populations[4]) -
           (populations[11] + populations[12]) -
           (populations[13] + populations[14]) -
           (populations[15] + populations[16]) -
           (populations[17] + populations[18]) -
           2.0f * ((populations[5] + populations[6]) -
                   (populations[7] + populations[8]) -
                   (populations[9] + populations[10]));
  case 7:
    return (populations[7] + populations[8]) -
           (populations[9] + populations[10]);
  case 8:
    return (populations[11] + populations[12]) -
           (populations[13] + populations[14]);
  case 9:
    return (populations[15] + populations[16]) -
           (populations[17] + populations[18]);
  case 10:
    return -2.0f * (populations[1] - populations[2]) +
           (populations[7] - populations[8]) +
           (populations[9] - populations[10]) +
           (populations[11] - populations[12]) +
           (populations[13] - populations[14]);
  case 11:
    return -2.0f * (populations[3] - populations[4]) +
           (populations[7] - populations[8]) -
           (populations[9] - populations[10]) +
           (populations[15] - populations[16]) +
           (populations[17] - populations[18]);
  case 12:
    return -2.0f * (populations[5] - populations[6]) +
           (populations[11] - populations[12]) -
           (populations[13] - populations[14]) +
           (populations[15] - populations[16]) -
           (populations[17] - populations[18]);
  case 13:
    return (populations[7] - populations[8]) +
           (populations[9] - populations[10]) -
           (populations[11] - populations[12]) -
           (populations[13] - populations[14]);
  case 14:
    return (populations[7] - populations[8]) -
           (populations[9] - populations[10]) -
           (populations[15] - populations[16]) -
           (populations[17] - populations[18]);
  case 15:
    return (populations[11] - populations[12]) -
           (populations[13] - populations[14]) -
           (populations[15] - populations[16]) +
           (populations[17] - populations[18]);
  case 16:
    return populations[0] + populations[7] + populations[8] + populations[9] +
           populations[10] + populations[11] + populations[12] +
           populations[13] + populations[14] + populations[15] +
           populations[16] + populations[17] + populations[18] -
           2.0f * ((populations[1] + populations[2]) +
                   (populations[3] + populations[4]) +
                   (populations[5] + populations[6]));
  case 17:
    return -(populations[1] + populations[2]) +
           (populations[3] + populations[4]) +
           (populations[11] + populations[12]) +
           (populations[13] + populations[14]) -
           (populations[15] + populations[16]) -
           (populations[17] + populations[18]);
  case 18:
    return -(populations[1] + populations[2]) -
           (populations[3] + populations[4]) -
           (populations[11] + populations[12]) -
           (populations[13] + populations[14]) -
           (populations[15] + populations[16]) -
           (populations[17] + populations[18]) +
           2.0f * ((populations[5] + populations[6]) +
                   (populations[7] + populations[8]) +
                   (populations[9] + populations[10]));
  }
  return 0.0;
}

/**
 *  @param[in]  node_index        Node index around (8) particle
 *  @param[out] mode              Local register values mode
 *  @param[in]  n_a               Local node residing in array a
 */
__device__ __inline__ void
calc_mass_and_momentum_mode(Utils::Array<float, 4> &mode, LB_nodes_gpu n_a,
                            unsigned int node_index) {
  /* mass mode */
  mode[0] = calc_mode_x_from_n(n_a.populations[node_index], 0);

  /* momentum modes */
  mode[1] = calc_mode_x_from_n(n_a.populations[node_index], 1);

  mode[2] = calc_mode_x_from_n(n_a.populations[node_index], 2);

  mode[3] = calc_mode_x_from_n(n_a.populations[node_index], 3);
}

struct LB_boundaries_gpu {
  /** For each fluid node this array contains either
   *  0 if the node is not a boundary, or the index of
   *  the boundary in LBBoundaries::lbboundaries minus one.
   */
  unsigned int *index = nullptr;
  /** If the node is a boundary node, this contains the
   *  velocity of the boundary
   */
  Utils::Array<float, 3> *velocity = nullptr;
};

inline __device__ float4 random_wrapper_philox(unsigned int index,
                                               unsigned int mode,
                                               uint64_t philox_counter) {
  // Split the 64 bit counter into two 32 bit ints.
  auto const philox_counter_hi = static_cast<uint32_t>(philox_counter >> 32);
  auto const philox_counter_low = static_cast<uint32_t>(philox_counter);
  uint4 rnd_ints =
      curand_Philox4x32_10(make_uint4(index, philox_counter_hi, 0, mode),
                           make_uint2(philox_counter_low, 0));
  float4 rnd_floats;
  rnd_floats.w = static_cast<float>(rnd_ints.w) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  rnd_floats.x = static_cast<float>(rnd_ints.x) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  rnd_floats.y = static_cast<float>(rnd_ints.y) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  rnd_floats.z = static_cast<float>(rnd_ints.z) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  return rnd_floats;
}

#endif // CUDA
#endif
