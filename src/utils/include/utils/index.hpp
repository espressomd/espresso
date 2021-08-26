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
#ifndef UTILS_INDEX_HPP
#define UTILS_INDEX_HPP

#include <cassert>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>

#include "Vector.hpp"

namespace Utils {

enum class MemoryOrder { COLUMN_MAJOR, ROW_MAJOR };

/** get the linear index from the position (@p a,@p b,@p c) in a 3D grid
 *  of dimensions @p adim.
 *
 * @return           The linear index
 * @param a , b , c  Position in 3D space
 * @param adim       Dimensions of the underlying grid
 * @param memory_order Row- or column-major
 */
inline int
get_linear_index(int a, int b, int c, const Vector3i &adim,
                 MemoryOrder memory_order = MemoryOrder::COLUMN_MAJOR) {
  assert((a >= 0) && (a < adim[0]));
  assert((b >= 0) && (b < adim[1]));
  assert((c >= 0) && (c < adim[2]));

  switch (memory_order) {
  case MemoryOrder::COLUMN_MAJOR:
    return a + adim[0] * (b + adim[1] * c);
  case MemoryOrder::ROW_MAJOR:
    return adim[1] * adim[2] * a + adim[2] * b + c;
  default:
    throw std::runtime_error("Unknown memory order");
  }
}

inline int
get_linear_index(const Vector3i &ind, const Vector3i &adim,
                 MemoryOrder memory_order = MemoryOrder::COLUMN_MAJOR) {
  return get_linear_index(ind[0], ind[1], ind[2], adim, memory_order);
}

/**
 * @brief Linear index into an upper triangular matrix.
 *
 * This is row-major.
 *
 * @tparam T Integral
 * @param i row index
 * @param j column index
 * @param n matrix size
 * @return linear index
 */
template <class T> T upper_triangular(T i, T j, T n) {
  /* n is a valid size */
  assert(n >= 0);
  /* i is a valid row index */
  assert((i >= 0) && (i < n));
  /* j is in the upper triangle */
  assert((j >= i) && (j < n));
  return (n * (n - 1)) / 2 - ((n - i) * (n - i - 1)) / 2 + j;
}

} // namespace Utils

#endif
