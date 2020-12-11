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

#include <boost/qvm/deduce_quat.hpp>
#include <boost/qvm/deduce_scalar.hpp>
#include <boost/qvm/deduce_vec.hpp>
#include <boost/qvm/quat.hpp>
#include <boost/qvm/quat_access.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/quat_traits.hpp>
#include <boost/qvm/quat_vec_operations.hpp>
#include <boost/serialization/array.hpp>

#include "utils/Array.hpp"
#include "utils/Vector.hpp"

namespace Utils {
template <typename T> struct Quaternion {
  using container = typename Utils::Array<T, 4>;
  container m_data;
  using pointer = typename container::pointer;
  using const_pointer = typename container::const_pointer;
  using value_type = typename container::value_type;
  using reference = typename container::reference;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &m_data;
  }

  void normalize() { boost::qvm::normalize(*this); }

  Quaternion<T> normalized() const { return boost::qvm::normalized(*this); }

  value_type operator[](std::size_t i) const { return m_data[i]; }
  reference operator[](std::size_t i) { return m_data[i]; }

  value_type norm() const { return boost::qvm::mag(*this); }
  value_type norm2() const { return boost::qvm::mag_sqr(*this); }

  static Quaternion<T> identity() { return boost::qvm::identity_quat<T>(); }

  static Quaternion<T> zero() { return boost::qvm::zero_quat<T>(); }
  constexpr pointer data() { return m_data.data(); }
  constexpr const_pointer data() const noexcept { return m_data.data(); }
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
  using scalar_type = T;

  template <std::size_t I>
  static constexpr inline scalar_type &write_element(quat_type &q) {
    static_assert(I < 4 and I >= 0, "Invalid index into quaternion.");
    return q[I];
  }

  template <std::size_t I>
  static constexpr inline scalar_type read_element(quat_type const &q) {
    static_assert(I < 4 and I >= 0, "Invalid index into quaternion.");
    return q[I];
  }

  static inline scalar_type read_element_idx(std::size_t i,
                                             quat_type const &q) {
    assert(i < 4 and i >= 0);
    return q[i];
  }

  static inline scalar_type &write_element_idx(std::size_t i, quat_type &q) {
    assert(i < 4 and i >= 0);
    return q[i];
  }
};

template <typename T, typename U>
struct deduce_scalar<Utils::Quaternion<T>, Utils::Quaternion<U>> {
  using type = std::common_type_t<T, U>;
};

template <typename T, typename U, std::size_t N>
struct deduce_scalar<Utils::Quaternion<T>, Utils::Vector<U, N>> {
  using type = std::common_type_t<T, U>;
};

template <typename T, typename U>
struct deduce_vec2<Utils::Quaternion<T>, Utils::Vector<U, 3>, 3> {
  using type = typename Utils::Vector<std::common_type_t<T, U>, 3>;
};

template <typename T> struct deduce_quat<Utils::Quaternion<T>> {
  using type = typename Utils::Quaternion<T>;
};

template <typename T, typename U>
struct deduce_quat2<Utils::Quaternion<T>, Utils::Quaternion<U>> {
  using type = typename Utils::Quaternion<std::common_type_t<T, U>>;
};

} // namespace qvm
} // namespace boost
#endif // SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
