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

#include <iterator>
#include <numeric>

#include "Vector.hpp"

namespace Utils {

/**
 * @brief Returns the flat index for given multidimensional indices and
 * dimensions.
 * @param unravelled_indices a container with the multidimensional
 * indices.
 * @param dimensions a container with the corresponding dimensions.
 * @retval the flat index
 */
template <typename T, typename U>
inline size_t ravel_index(const T &unravelled_indices, const U &dimensions) {
  const auto n_dims = unravelled_indices.size();
  if (n_dims != dimensions.size()) {
    throw std::invalid_argument(
        "Index vector and dimensions vector must have same dimensions.");
  }
  std::size_t res = unravelled_indices.back();
  std::size_t temp_prod = 1;
  for (int i = unravelled_indices.size() - 2; i >= 0; --i) {
    temp_prod *= dimensions[i + 1];
    res += unravelled_indices[i] * temp_prod;
  }
  return res;
}

/**
 * @brief Returns the unraveled index of the provided flat index.
 *        Therefore is the inversion of flattening an ndims dimensional index.
 * @param[in] dimensions_begin Iterator pointing to the begin of the container
 * with the lengths of the dimensions.
 * @param[in] dimensions_end   Iterator pointing to the end of the container
 * with the lengths of the dimensions.
 * @param[out] begin_out       Outputiterator pointing to the begin of the
 * container where the result should be written to.
 * @param[in] ravelled_index   The flat index.
 */
template <typename InputIterator, typename OutputIterator, typename T>
inline void unravel_index(InputIterator dimensions_begin,
                          InputIterator dimensions_end,
                          OutputIterator begin_out, T ravelled_index) {
  auto end_out = begin_out + std::distance(dimensions_begin, dimensions_end);
  auto rbegin_in = std::make_reverse_iterator(dimensions_begin);
  auto rend_in = std::make_reverse_iterator(dimensions_end);
  auto rend_out = std::make_reverse_iterator(end_out);
  std::size_t mul = 1;
  for (; rend_in != rbegin_in; ++rend_in, ++rend_out) {
    *rend_out = (ravelled_index / mul) % (*rend_in);
    mul *= (*rend_in);
  }
}

enum class MemoryOrder { COLUMN_MAJOR, ROW_MAJOR };

/*************************************************************/
/** \name Three dimensional grid operations                  */
/*************************************************************/
/*@{*/

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

/*@}*/

} // namespace Utils

#endif
