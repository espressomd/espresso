/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <cassert>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>

#include "Vector.hpp"
#include "device_qualifier.hpp"

namespace Utils {

enum class MemoryOrder { COLUMN_MAJOR, ROW_MAJOR };

template <MemoryOrder memory_order>
DEVICE_QUALIFIER inline int get_linear_index(int a, int b, int c,
                                             Vector3i const &adim) {
  DEVICE_ASSERT((a >= 0) && (a < adim[0]));
  DEVICE_ASSERT((b >= 0) && (b < adim[1]));
  DEVICE_ASSERT((c >= 0) && (c < adim[2]));
  if constexpr (memory_order == MemoryOrder::COLUMN_MAJOR) {
    return a + adim[0] * (b + adim[1] * c);
  }
  return adim[1] * adim[2] * a + adim[2] * b + c;
}

/** Get the linear index from the position (@p a,@p b,@p c) in a 3D grid
 *  of dimensions @p adim.
 *
 * @return           The linear index
 * @param a , b , c  Position in 3D space
 * @param adim       Dimensions of the underlying grid
 * @param memory_order Row- or column-major
 */
DEVICE_QUALIFIER inline int
get_linear_index(int a, int b, int c, Vector3i const &adim,
                 MemoryOrder memory_order = MemoryOrder::COLUMN_MAJOR) {
  if (memory_order == MemoryOrder::COLUMN_MAJOR) {
    return get_linear_index<MemoryOrder::COLUMN_MAJOR>(a, b, c, adim);
  }
  return get_linear_index<MemoryOrder::ROW_MAJOR>(a, b, c, adim);
}

DEVICE_QUALIFIER inline int
get_linear_index(Vector3i const &ind, Vector3i const &adim,
                 MemoryOrder memory_order = MemoryOrder::COLUMN_MAJOR) {
  return get_linear_index(ind[0], ind[1], ind[2], adim, memory_order);
}

template <MemoryOrder memory_order>
DEVICE_QUALIFIER inline int get_linear_index(Vector3i const &ind,
                                             Vector3i const &adim) {
  return get_linear_index<memory_order>(ind[0], ind[1], ind[2], adim);
}

/**
 * @brief Linear index into a lower triangular matrix.
 *
 * This is row-major.
 *
 * @tparam T Integral type
 * @param i row index
 * @param j column index
 * @return linear index
 */
template <class T> DEVICE_QUALIFIER T lower_triangular(T i, T j) {
  /* i is a valid row index */
  DEVICE_ASSERT(i >= 0);
  /* j is in the lower triangle */
  DEVICE_ASSERT(j >= 0 and j <= i);
  return (i * (i + 1)) / 2 + j;
}

} // namespace Utils
