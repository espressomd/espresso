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
 *  Implementation in lbgpu_cuda.cu.
 */

#ifndef LBGPU_CUH
#define LBGPU_CUH

#include "config.hpp"

#ifdef CUDA
#include "curand_wrapper.hpp"

#include <utils/Array.hpp>

#ifdef CUDA
/** Velocity densities for the lattice Boltzmann system. */
struct LB_nodes_gpu {
  /** velocity density of the node */
  float *vd = nullptr;
  unsigned int *boundary = nullptr;
  Utils::Array<float, 3> *boundary_velocity = nullptr;
};

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
#endif //  CUDA

#endif // CUDA
#endif
