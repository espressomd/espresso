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
#ifndef UTILS_MATH_TENSOR_PRODUCT_HPP
#define UTILS_MATH_TENSOR_PRODUCT_HPP

#include "utils/Vector.hpp"
#include "utils/matrix.hpp"

#include <algorithm>
#include <cstddef>

namespace Utils {
template <typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> tensor_product(const Vector<T, N> &x, const Vector<T, M> &y) {
  return boost::qvm::col_mat(x) * boost::qvm::row_mat(y);
}

/*
 * @brief Overload if left operand is scalar.
 */
template <typename T, std::size_t N>
Vector<T, N> tensor_product(const T &x, const Vector<T, N> &y) {
  return x * y;
}
} // namespace Utils

#endif
