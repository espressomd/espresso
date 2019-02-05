/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef CORE_EXTERNAL_FIELD_FIELDS_AFFINE_MAP_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_AFFINE_MAP_HPP

#include "jacobian_type.hpp"
#include "utils/Vector.hpp"

namespace FieldCoupling {
namespace Fields {

namespace detail {

template <typename T, size_t codim> struct matrix_vector_impl {
  Vector<codim, T> operator()(const Vector<codim, Vector<3, T>> &A,
                              Vector<3, T> const &v) const {
    Vector<codim, T> ret;

    for (int i = 0; i < codim; i++)
      ret[i] = A[i] * v;

    return ret;
  }
};

template <typename T> struct matrix_vector_impl<T, 1> {
  T operator()(const Vector<3, T> &A, Vector<3, T> const &v) const {
    return A * v;
  }
};
} // namespace detail

/**
 * @brief Affine transform of an vector field.
 *
 * Returns A * x + b, where * is matrix multiplication.
 */
template <typename T, size_t codim> class AffineMap {
public:
  using value_type = typename decay_to_scalar<Vector<codim, T>>::type;
  using jacobian_type = detail::jacobian_type<T, codim>;

private:
  jacobian_type m_A;
  value_type m_b;

public:
  AffineMap(const jacobian_type &A, const value_type &b) : m_A(A), m_b(b) {}

  jacobian_type &A() { return m_A; }
  value_type &b() { return m_b; }

  value_type operator()(const Vector3d &pos, double = {}) const {
    return detail::matrix_vector_impl<T, codim>{}(m_A, pos) + m_b;
  }

  jacobian_type jacobian(const Vector3d &, double = {}) const { return m_A; }

  bool fits_in_box(const Vector3d &) const { return true; }
};
} // namespace Fields
} // namespace FieldCoupling

#endif
