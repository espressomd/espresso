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
#ifndef CURAND_WRAPPER_HPP
#define CURAND_WRAPPER_HPP

#if defined(__CUDACC__)

#include <curand_kernel.h>

#elif defined(__HIPCC__)

#include <rocrand/rocrand_kernel.h>

#define CURAND_2POW32_INV ROCRAND_2POW32_INV

class philox4x32_10_stateless : private rocrand_device::philox4x32_10_engine {
public:
  FQUALIFIERS
  philox4x32_10_stateless() {}

  FQUALIFIERS
  uint4 operator()(uint4 counter, uint2 key) {
    return ten_rounds(counter, key);
  }
};

__forceinline__ __device__ uint4 curand_Philox4x32_10(uint4 counter,
                                                      uint2 key) {
  philox4x32_10_stateless e;
  return e(counter, key);
}

#endif

#endif // CURAND_WRAPPER_HPP
