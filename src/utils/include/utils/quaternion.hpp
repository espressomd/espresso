/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
#define SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP

#include <boost/qvm/deduce_scalar.hpp>
#include <boost/qvm/deduce_vec.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/quat_traits.hpp>
#include <boost/qvm/quat_vec_operations.hpp>

#include "utils/Vector.hpp"

namespace Utils {
template <typename T> class Quaternion : private Vector<T, 4> {
  using base = Vector<T, 4>;
  friend class boost::serialization::access;

public:
  using base::base;
  using base::operator[];
  using base::begin;
  using base::data;
  using base::end;
  using typename base::value_type;
  Quaternion<T>() = default;
  Quaternion<T>(Quaternion<T> const &) = default;
  Quaternion<T> &operator=(Quaternion<T> const &) = default;
  Quaternion<T>(std::initializer_list<T> values) {
    std::copy(values.begin(), values.end(), begin());
  };

  Quaternion &normalize() {
    boost::qvm::normalize(*this);
    return *this;
  }

  T norm() const { return boost::qvm::mag(*this); }
  T norm2() const { return boost::qvm::mag_sqr(*this); }

  static Quaternion<T> identity() { return boost::qvm::identity_quat<T>(); }

  static Quaternion<T> zero() { return boost::qvm::zero_quat<T>(); }
};

template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic<U>::value, bool> = true>
Quaternion<T> operator*(const U &b, const Quaternion<T> &a) {
  return boost::qvm::operator*(a, b);
}

using boost::qvm::operator*;
using boost::qvm::operator*=;
using boost::qvm::operator==;
using boost::qvm::dot;
using boost::qvm::operator/;
using boost::qvm::operator/=;
using boost::qvm::operator+;
using boost::qvm::operator+=;
using boost::qvm::operator-;
using boost::qvm::operator-=;

} // namespace Utils

namespace boost {

namespace qvm {

template <class T> struct quat_traits<Utils::Quaternion<T>> {
  using quat_type = Utils::Quaternion<T>;
  using scalar_type = typename quat_type::value_type;

  template <std::size_t I>
  static constexpr inline scalar_type &write_element(quat_type &q) {
    return q[I];
  }

  template <std::size_t I>
  static constexpr inline scalar_type read_element(quat_type const &q) {
    return q[I];
  }
};

template <typename T, typename U>
struct deduce_scalar<Utils::Quaternion<T>, Utils::Quaternion<U>> {
  using type = std::common_type_t<T, U>;
};

template <typename T, typename U>
struct deduce_vec2<Utils::Quaternion<T>, Utils::Vector<U, 3>, 3> {
  using type = typename Utils::Vector<std::common_type_t<T, U>, 3>;
};

} // namespace qvm
} // namespace boost
#endif // SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
