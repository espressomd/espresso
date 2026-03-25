/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#pragma once

#include <gpu/FieldAccessor.h>

#if defined(__CUDACC__)
#include <cufft.h>
#endif

#include <cmath>
#include <numbers>
#include <type_traits>

namespace walberla {

template <typename FloatType>
FloatType greens_function(int x, int y, int z, auto const &dim) {
  if (x == 0 && y == 0 && z == 0)
    return 0.;
  auto constexpr two_pi = FloatType{2} * std::numbers::pi_v<FloatType>;
  return FloatType(-0.5) /
         FloatType(std::cos(two_pi * FloatType(x) / FloatType(dim[0])) +
                   std::cos(two_pi * FloatType(y) / FloatType(dim[1])) +
                   std::cos(two_pi * FloatType(z) / FloatType(dim[2])) - 3.) /
         static_cast<FloatType>(dim[0] * dim[1] * dim[2]);
}

#if defined(__CUDACC__)
// LCOV_EXCL_START
template <typename FloatType>
__global__ void
create_greens_function(gpu::FieldAccessor<FloatType> greens_function, int x_min,
                       int y_min, int z_min, int x_max, int y_max, int z_max,
                       int global_dim_x, int global_dim_y, int global_dim_z) {
  using RealType = std::conditional<std::is_same<FloatType, float>::value,
                                    cufftReal, cufftDoubleReal>::type;
  greens_function.set(blockIdx, threadIdx);
  unsigned int index =
      greens_function.getLinearIndex(blockIdx, threadIdx, gridDim, blockDim);
  unsigned int local_dim[3] = {static_cast<unsigned int>(x_max - x_min),
                               static_cast<unsigned int>(y_max - y_min),
                               static_cast<unsigned int>(z_max - z_min)};
  unsigned int yz_slice_size = index / local_dim[0];
  auto global_coord_x = x_min + index % local_dim[0];
  auto global_coord_y = y_min + yz_slice_size % local_dim[1];
  auto global_coord_z = z_min + yz_slice_size / local_dim[1];
  constexpr RealType two_pi{2. * M_PI};
  if (index < local_dim[0] * local_dim[1] * local_dim[2]) {
    if (index == 0u and x_min == 0 and y_min == 0 and z_min == 0) {
      // setting 0th Fourier mode to 0 enforces charge neutrality
      greens_function.get(0u) = RealType{0};
    } else {
      greens_function.get(0u) =
          RealType{-0.5} /
          (cos(two_pi * static_cast<RealType>(global_coord_x) /
               static_cast<RealType>(global_dim_x)) +
           cos(two_pi * static_cast<RealType>(global_coord_y) /
               static_cast<RealType>(global_dim_y)) +
           cos(two_pi * static_cast<RealType>(global_coord_z) /
               static_cast<RealType>(global_dim_z)) -
           RealType{3}) /
          static_cast<RealType>(global_dim_x * global_dim_y * global_dim_z);
    }
  }
}

template <typename FloatType, typename ComplexType>
__global__ void
multiply_by_greens_function(gpu::FieldAccessor<ComplexType> potential,
                            gpu::FieldAccessor<FloatType> greens_function) {
  potential.set(blockIdx, threadIdx);
  greens_function.set(blockIdx, threadIdx);
  if (potential.isValidPosition() && greens_function.isValidPosition()) {
    potential.get(0u) =
        ComplexType(potential.get(0u).x * greens_function.get(0u),
                    potential.get(0u).y * greens_function.get(0u));
  }
}

template <typename FloatType>
__global__ void add_fields_with_factor(gpu::FieldAccessor<FloatType> field_out,
                                       gpu::FieldAccessor<FloatType> field_add,
                                       FloatType const factor) {
  field_out.set(blockIdx, threadIdx);
  field_add.set(blockIdx, threadIdx);
  if (field_out.isValidPosition() && field_add.isValidPosition()) {
    field_out.get(0u) += field_add.get(0u) * factor;
  }
}

template <typename FloatType>
__global__ void move_field(gpu::FieldAccessor<FloatType> dest_field,
                           gpu::FieldAccessor<FloatType> src_field) {
  dest_field.set(blockIdx, threadIdx);
  src_field.set(blockIdx, threadIdx);
  if (dest_field.isValidPosition() && src_field.isValidPosition()) {
    dest_field.get(0u) = src_field.get(0u);
    src_field.get(0u) = FloatType{0};
  }
}
// LCOV_EXCL_STOP
#endif // __CUDACC__

} // namespace walberla
