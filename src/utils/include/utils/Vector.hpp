/*
 * Copyright (C) 2014-2019 The ESPResSo project
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

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <vector>

#include "utils/Array.hpp"

namespace Utils {

template <typename T, std::size_t N> class Vector : public Array<T, N> {
  using Base = Array<T, N>;

public:
  using Base::Base;
  using Array<T, N>::at;
  using Array<T, N>::operator[];
  using Array<T, N>::front;
  using Array<T, N>::back;
  using Array<T, N>::data;
  using Array<T, N>::begin;
  using Array<T, N>::cbegin;
  using Array<T, N>::end;
  using Array<T, N>::cend;
  using Array<T, N>::empty;
  using Array<T, N>::size;
  using Array<T, N>::max_size;
  using Array<T, N>::fill;
  using Array<T, N>::broadcast;
  Vector() = default;
  Vector(Vector const &) = default;
  Vector &operator=(Vector const &) = default;

  void swap(Vector &rhs) { std::swap_ranges(begin(), end(), rhs.begin()); }

private:
  constexpr void copy_init(T const *first, T const *last) {
    auto it = begin();
    while (first != last) {
      *it++ = *first++;
    }
  }

public:
  template <class Range>
  explicit Vector(Range const &rng) : Vector(std::begin(rng), std::end(rng)) {}
  explicit constexpr Vector(T const (&v)[N]) : Base() {
    copy_init(std::begin(v), std::end(v));
  }

  constexpr Vector(std::initializer_list<T> v) : Base() {
    if (N != v.size()) {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }

    copy_init(v.begin(), v.end());
  }

  template <typename InputIterator>
  Vector(InputIterator first, InputIterator last) : Base() {
    if (std::distance(first, last) == N) {
      std::copy_n(first, N, begin());
    } else {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }
  }

  /**
   * @brief Create a vector that has all entries set to
   *         one value.
   */
  static Vector<T, N> broadcast(T const &s) {
    Vector<T, N> ret;
    std::fill(ret.begin(), ret.end(), s);

    return ret;
  }

  std::vector<T> as_vector() const { return std::vector<T>(begin(), end()); }

  operator std::vector<T>() const { return as_vector(); }

  template <class U> explicit operator Vector<U, N>() const {
    Vector<U, N> ret;

    std::transform(begin(), end(), ret.begin(),
                   [](auto e) { return static_cast<U>(e); });

    return ret;
  }

  inline T norm2() const { return (*this) * (*this); }
  inline T norm() const { return std::sqrt(norm2()); }

  /*
   * @brief Normalize the vector.
   *
   * Normalize the vector by its length,
   * if not zero, otherwise the vector is unchanged.
   */

  inline Vector &normalize() {
    auto const l = norm();
    if (l > T(0)) {
      for (int i = 0; i < N; i++)
        this->operator[](i) /= l;
    }

    return *this;
  }
};

template <class T> using Vector3 = Vector<T, 3>;

template <size_t N> using VectorXd = Vector<double, N>;
using Vector2d = VectorXd<2>;
using Vector3d = VectorXd<3>;
using Vector4d = VectorXd<4>;
using Vector6d = VectorXd<6>;
using Vector9d = VectorXd<9>;
using Vector19d = VectorXd<19>;

template <size_t N> using VectorXf = Vector<float, N>;
using Vector3f = VectorXf<3>;

template <size_t N> using VectorXi = Vector<int, N>;
using Vector3i = VectorXi<3>;

template <class T, size_t N, size_t M> using Matrix = Vector<Vector<T, M>, N>;

/**
 * @brief Trace of a matrix.
 *
 * Returns the sum of the diagonal elements
 * of a square matrix.
 *
 * @tparam T Arithmetic type
 * @tparam N Matrix dimension
 * @param m Input matrix
 * @return Trace of matrix.
 */
template <class T, size_t N> T trace(Matrix<T, N, N> const &m) {
  auto tr = T{};
  for (size_t i = 0; i < N; i++)
    tr += m[i][i];

  return tr;
}

/**
 * @brief Flatten a matrix to a linear vector.
 *
 * @param m Input Matrix
 * @return Flat vector with elements of the matrix.
 */
template <class T, size_t N, size_t M>
Vector<T, N * M> flatten(Matrix<T, N, M> const &m) {
  Vector<T, N * M> ret;

  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < M; j++)
      ret[i * M + j] = m[j][i];

  return ret;
}

namespace detail {
template <size_t N, typename T, typename U, typename Op>
auto binary_op(Vector<T, N> const &a, Vector<U, N> const &b, Op op) {
  using std::declval;

  using R = decltype(op(declval<T>(), declval<U>()));
  Vector<R, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(ret),
                 op);

  return ret;
}

template <size_t N, typename T, typename Op>
Vector<T, N> &binary_op_assign(Vector<T, N> &a, Vector<T, N> const &b, Op op) {
  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(a), op);
  return a;
}

template <size_t N, typename T, typename Op>
constexpr bool all_of(Vector<T, N> const &a, Vector<T, N> const &b, Op op) {
  for (int i = 0; i < a.size(); i++) {
    /* Short circuit */
    if (!static_cast<bool>(op(a[i], b[i]))) {
      return false;
    }
  }

  return true;
}
} // namespace detail

template <size_t N, typename T>
constexpr bool operator<(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::less<T>());
}

template <size_t N, typename T>
constexpr bool operator>(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::greater<T>());
}

template <size_t N, typename T>
constexpr bool operator<=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::less_equal<T>());
}

template <size_t N, typename T>
constexpr bool operator>=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::greater_equal<T>());
}

template <size_t N, typename T>
constexpr bool operator==(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::equal_to<T>());
}

template <size_t N, typename T>
constexpr bool operator!=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return not(a == b);
}

template <size_t N, typename T, typename U>
auto operator+(Vector<T, N> const &a, Vector<U, N> const &b) {
  return detail::binary_op(a, b, std::plus<>());
}

template <size_t N, typename T>
Vector<T, N> &operator+=(Vector<T, N> &a, Vector<T, N> const &b) {
  return detail::binary_op_assign(a, b, std::plus<T>());
}

template <size_t N, typename T, typename U>
auto operator-(Vector<T, N> const &a, Vector<U, N> const &b) {
  return detail::binary_op(a, b, std::minus<>());
}

template <size_t N, typename T> Vector<T, N> operator-(Vector<T, N> const &a) {
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(ret),
                 [](T const &v) { return -v; });

  return ret;
}

template <size_t N, typename T>
Vector<T, N> &operator-=(Vector<T, N> &a, Vector<T, N> const &b) {
  return detail::binary_op_assign(a, b, std::minus<T>());
}

/* Scalar multiplication */
template <size_t N, typename T, class U>
auto operator*(U const &a, Vector<T, N> const &b) {
  using R = decltype(a * std::declval<T>());
  Vector<R, N> ret;

  std::transform(std::begin(b), std::end(b), std::begin(ret),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T, class U>
auto operator*(Vector<T, N> const &b, U const &a) {
  using R = decltype(std::declval<T>() * a);
  Vector<R, N> ret;

  std::transform(std::begin(b), std::end(b), std::begin(ret),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T>
Vector<T, N> &operator*=(Vector<T, N> &b, T const &a) {
  std::transform(std::begin(b), std::end(b), std::begin(b),
                 [a](T const &val) { return a * val; });
  return b;
}

/* Scalar division */
template <size_t N, typename T>
Vector<T, N> operator/(Vector<T, N> const &a, T const &b) {
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), ret.begin(),
                 [b](T const &val) { return val / b; });
  return ret;
}

template <size_t N, typename T>
Vector<T, N> &operator/=(Vector<T, N> &a, T const &b) {
  std::transform(std::begin(a), std::end(a), std::begin(a),
                 [b](T const &val) { return val / b; });
  return a;
}

namespace detail {
template <class T> struct is_vector : std::false_type {};
template <class T, size_t N> struct is_vector<Vector<T, N>> : std::true_type {};
} // namespace detail

/* Scalar product */
template <size_t N, typename T, class U,
          class = std::enable_if_t<not(detail::is_vector<T>::value or
                                       detail::is_vector<U>::value)>>
auto operator*(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  return std::inner_product(std::begin(a), std::end(a), std::begin(b), R{});
}

/* Matrix x Vector */
template <class T>
Vector<T, 3> operator*(Matrix<T, 3, 3> const &A, Vector<T, 3> const &v) {
  return {v[0] * A[0][0] + v[1] * A[0][1] + v[2] * A[0][2],
          v[0] * A[1][0] + v[1] * A[1][1] + v[2] * A[1][2],
          v[0] * A[2][0] + v[1] * A[2][1] + v[2] * A[2][2]};
}

/**
 * @brief Transpose matrix
 * @param m Input matrix
 * @return Transposed matrix
 */
template <class T> Matrix<T, 3, 3> transpose(Matrix<T, 3, 3> const &m) {
  return {
      {m[0][0], m[1][0], m[2][0]},
      {m[0][1], m[1][1], m[2][1]},
      {m[0][2], m[1][2], m[2][2]},
  };
}

/**
 * @brief Diagonal matrix with diagonal elements from vector.
 *
 * Diagonal matrix with vector entries as diagonal:
 *
 * \f[
 *     D_{ij} = \delta_{ij} v_i
 * \f]
 *
 * Only implemented for 3x3 matrices.
 *
 * @tparam T scalar type
 * @param v Vector with diagonal elements
 */
template <class T> Matrix<T, 3, 3> diag_matrix(Vector<T, 3> const &v) {
  return {{v[0], 0, 0}, {0, v[1], 0}, {0, 0, v[2]}};
}

/**
 * @brief Matrix product
 *
 * Matrix product C, where
 *
 * \f[
 *     C_{ij} = \sum_k A_{ik} B_{kj}
 * \f]
 *
 * Only implemented for 3x3 matrices.
 *
 * @tparam T scalar type
 * @param A Left-hand side
 * @param B Right-hand side
 * @return Matrix product
 */
template <class T>
Matrix<T, 3, 3> operator*(Matrix<T, 3, 3> const &A, Matrix<T, 3, 3> const &B) {
  auto const Bt = transpose(B);

  return {{A[0] * Bt[0], A[0] * Bt[1], A[0] * Bt[2]},
          {A[1] * Bt[0], A[1] * Bt[1], A[1] * Bt[2]},
          {A[2] * Bt[0], A[2] * Bt[1], A[2] * Bt[2]}};
}

template <size_t N, typename T, class U,
          class = std::enable_if_t<std::is_integral<T>::value &&
                                   std::is_integral<U>::value>>
auto operator%(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() % declval<U>());
  Vector<R, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(ret),
                 [](T const &ai, U const &bi) { return ai % bi; });

  return ret;
}

/* Componentwise square root */
template <size_t N, typename T> Vector<T, N> sqrt(Vector<T, N> const &a) {
  using std::sqrt;
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), ret.begin(),
                 [](T const &v) { return sqrt(v); });

  return ret;
}

template <class T>
Vector<T, 3> vector_product(Vector<T, 3> const &a, Vector<T, 3> const &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

template <class T, class U, size_t N>
auto hadamard_product(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret;
  std::transform(a.cbegin(), a.cend(), b.cbegin(), ret.begin(),
                 [](auto ai, auto bi) { return ai * bi; });

  return ret;
}

// specializations for when one or both operands is a scalar depending on
// compile time features (e.g. when PARTICLE_ANISOTROPY is not enabled)
template <class T, class U, size_t N,
          class = std::enable_if_t<not(detail::is_vector<T>::value)>>
auto hadamard_product(T const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret = a * b;

  return ret;
}

template <class T, class U, size_t N,
          class = std::enable_if_t<not(detail::is_vector<U>::value)>>
auto hadamard_product(Vector<T, N> const &a, U const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret = a * b;

  return ret;
}

template <typename T, typename U,
          class = std::enable_if_t<not(detail::is_vector<T>::value or
                                       detail::is_vector<U>::value)>>
auto hadamard_product(T const &a, U const &b) {
  return a * b;
}

template <class T, class U, size_t N>
auto hadamard_division(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret;
  std::transform(a.cbegin(), a.cend(), b.cbegin(), ret.begin(),
                 [](auto ai, auto bi) { return ai / bi; });

  return ret;
}

// specializations for when one or both operands is a scalar depending on
// compile time features (e.g. when PARTICLE_ANISOTROPY is not enabled)
template <class T, class U, size_t N,
          class = std::enable_if_t<not(detail::is_vector<U>::value)>>
auto hadamard_division(Vector<T, N> const &a, U const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret = a / b;

  return ret;
}

template <class T, class U, size_t N,
          class = std::enable_if_t<not(detail::is_vector<T>::value)>>
auto hadamard_division(T const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  Vector<R, N> ret;
  std::transform(std::begin(b), std::end(b), ret.begin(),
                 [a](T const &bi) { return a / bi; });
  return ret;
}

template <typename T, typename U,
          class = std::enable_if_t<not(detail::is_vector<T>::value or
                                       detail::is_vector<U>::value)>>
auto hadamard_division(T const &a, U const &b) {
  return a / b;
}

/**
 * @brief Meta function to turns a Vector<1, T> into T.
 */
template <typename T> struct decay_to_scalar {};
template <typename T, size_t N> struct decay_to_scalar<Vector<T, N>> {
  using type = Vector<T, N>;
};

template <typename T> struct decay_to_scalar<Vector<T, 1>> { using type = T; };

template <std::size_t I, class T, std::size_t N>
struct tuple_element<I, Vector<T, N>> {
  using type = T;
};

template <class T, std::size_t N>
struct tuple_size<Vector<T, N>> : std::integral_constant<std::size_t, N> {};

template <std::size_t I, class T, std::size_t N>
auto get(Vector<T, N> const &a) -> std::enable_if_t<(I < N), const T &> {
  return a[I];
}
} // namespace Utils
#endif
