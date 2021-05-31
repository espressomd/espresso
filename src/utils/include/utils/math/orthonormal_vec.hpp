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
#ifndef ESPRESSO_ORTHONORMAL_VEC_HPP
#define ESPRESSO_ORTHONORMAL_VEC_HPP

#include "utils/Vector.hpp"
#include "utils/constants.hpp"

namespace Utils {
/**
 * @brief Return a vector that is orthonormal to vec
 */
template <typename T, std::size_t N>
Vector<T, N> calc_orthonormal_vector(Vector<T, N> const &vec) {
  /* Calculate orthonormal vector using Gram-Schmidt orthogonalization of a
   trial vector. Only works if the trial vector is not parallel, so we have to
   try a second one in that case
  */
  Vector<Vector<T, N>, 2> try_vectors = {Vector<T, N>::broadcast(0),
                                         Vector<T, N>::broadcast(0)};
  try_vectors[0][0] = 1;
  try_vectors[1][1] = 1;

  Vector<T, N> ret;
  for (auto v : try_vectors) {
    auto orth_component = v - (v * vec) / vec.norm2() * vec;
    auto norm = orth_component.norm();
    if (norm >= 1. / Utils::sqrt_2()) {
      ret = orth_component / norm;
      break;
    }
  }
  return ret;
}

} // namespace Utils

#endif // ESPRESSO_ORTHONORMAL_VEC_HPP