/*
  Copyright (C) 2016,2017 The ESPResSo project

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
#ifndef CORE_UTILS_HISTOGRAM_HPP
#define CORE_UTILS_HISTOGRAM_HPP
namespace Utils {

/**
 * \brief Returns the unravelled index of the provided flat index.
 *        Therefore is the inversion of flattening an ndims dimensional index.
 * @param len_dims an int array of length ndims containing the lengths of the
 * dimensions. (Input)
 * @param ndims int denoting the number of dimensions. (Input)
 * @flattened_index an int denoting the flat index. (Input)
 * @unravelled_index_out an int array with length ndims where the unflat indices
 * are written to. (Output)
 */
inline void unravel_index(const int *const len_dims, const int ndims,
                          const int flattened_index,
                          int *unravelled_index_out) {
  // idea taken from
  // http://codinghighway.com/2014/02/22/c-multi-dimensional-arrays-part-2-flattened-to-unflattened-index/
  std::vector<int> mul(ndims);
  mul[ndims - 1] = 1;
  for (int j = ndims - 2; j >= 0; j--)
    mul[j] = mul[j + 1] * len_dims[j + 1];
  for (int j = 0; j < ndims; j++)
    unravelled_index_out[j] = (flattened_index / mul[j]) % len_dims[j];
}
} // Namespace Utils
#endif // CORE_UTILS_HISTOGRAM_HPP
