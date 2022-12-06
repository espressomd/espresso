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
#ifndef SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
#define SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP

/**
 * @file
 *
 * @brief Quaternion implementation and trait types
 * for boost qvm interoperability.
 */

#include <boost/qvm/deduce_quat.hpp>
#include <boost/qvm/deduce_scalar.hpp>
#include <boost/qvm/deduce_vec.hpp>
#include <boost/qvm/quat.hpp>
#include <boost/qvm/quat_access.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/quat_traits.hpp>
#include <boost/qvm/quat_vec_operations.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/level.hpp>

#include "utils/Array.hpp"
#include "utils/Vector.hpp"
#include "utils/matrix.hpp"

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace Utils {

/**
 * Quaternion representation.
 * @tparam T Element data type.
 */
template <typename T> struct Quaternion {
  using container = typename Utils::Array<T, 4>;
  container m_data;
  using pointer = typename container::pointer;
  using const_pointer = typename container::const_pointer;
  using value_type = typename container::value_type;
  using reference = typename container::reference;

private:
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &m_data;
  }

public:
  /**
   * @brief Normalize the quaternion in place.
   */
  void normalize() { boost::qvm::normalize(*this); }

  /**
   * @brief Retrieve a normalized copy of the quaternion.
   * @return Normalized quaternion.
   */
  Quaternion<T> normalized() const { return boost::qvm::normalized(*this); }

  /**
   * @brief Element access (const).
   * @param i Element index.
   * @return Value of element @p i.
   */
  value_type operator[](std::size_t i) const { return m_data[i]; }
  /**
   * @brief Element access (non const).
   * @param i Element index.
   * @return Value of element @p i.
   */
  reference operator[](std::size_t i) { return m_data[i]; }

  /**
   * @brief Retrieve the norm of the quaternion.
   * @return The norm.
   */
  value_type norm() const { return boost::qvm::mag(*this); }
  /**
   * @brief Retrieve the square of the norm of the quaternion.
   * @return The squared norm.
   */
  value_type norm2() const { return boost::qvm::mag_sqr(*this); }

  /**
   * @brief Construct an identity quaternion.
   * @return Identity quaternion.
   */
  static Quaternion<T> identity() { return boost::qvm::identity_quat<T>(); }

  /**
   * @brief Construct a zero quaternion.
   * @return Quaternion with all elements set to zero.
   */
  static Quaternion<T> zero() { return boost::qvm::zero_quat<T>(); }

  /**
   * @brief Access to the underlying data (non const).
   * @return Pointer to the first data member.
   */
  constexpr pointer data() { return m_data.data(); }

  /**
   * @brief Access to the underlying data (const).
   * @return Pointer to the first data member.
   */
  constexpr const_pointer data() const noexcept { return m_data.data(); }
};

/**
 * @brief Product quaternion and arithmetic type.
 * @tparam T Data type of quaternion @p a.
 * @tparam U Type of multiplier @p b.
 * @param b Quaternion.
 * @param a Multiplier.
 * @return Multiplied quaternion.
 */
template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
Quaternion<T> operator*(const U &b, const Quaternion<T> &a) {
  return boost::qvm::operator*(a, b);
}

/**
 * @brief Convert quaternion to rotation matrix.
 * @tparam T Data type of quaternion.
 * @param q Quaternion.
 * @return Rotation matrix.
 */
template <typename T> Matrix<T, 3, 3> rotation_matrix(Quaternion<T> const &q) {
  auto const normed_q = q.normalized();
  auto const id_mat = Utils::identity_mat<double, 3, 3>();
  auto const v1 = normed_q * id_mat.col<0>();
  auto const v2 = normed_q * id_mat.col<1>();
  auto const v3 = normed_q * id_mat.col<2>();
  return {{v1[0], v2[0], v3[0]}, {v1[1], v2[1], v3[1]}, {v1[2], v2[2], v3[2]}};
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
    static_assert(I < 4, "Invalid index into quaternion.");
    return q[I];
  }

  template <std::size_t I>
  static constexpr inline scalar_type read_element(quat_type const &q) {
    static_assert(I < 4, "Invalid index into quaternion.");
    return q[I];
  }

  static inline scalar_type read_element_idx(std::size_t i,
                                             quat_type const &q) {
    assert(i < 4);
    return q[i];
  }

  static inline scalar_type &write_element_idx(std::size_t i, quat_type &q) {
    assert(i < 4);
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

UTILS_ARRAY_BOOST_MPI_T(Utils::Quaternion, 0)
UTILS_ARRAY_BOOST_BIT_S(Utils::Quaternion, 0)
UTILS_ARRAY_BOOST_CLASS(Utils::Quaternion, 0, object_serializable)
UTILS_ARRAY_BOOST_TRACK(Utils::Quaternion, 0, track_never)

#endif // SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
