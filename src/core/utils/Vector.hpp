/*
  Copyright (C) 2014-2018 The ESPResSo project

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
  Vector(InputIterator first, InputIterator last) {
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

template <size_t N> using VectorXd = Vector<double, N>;
using Vector2d = VectorXd<2>;
using Vector3d = VectorXd<3>;
using Vector4d = VectorXd<4>;
using Vector6d = VectorXd<6>;
using Vector19d = VectorXd<19>;

using Vector3i = Vector<int, 3>;

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

/* Scalar product */
template <size_t N, typename T, class U>
auto operator*(Vector<T, N> const &a, Vector<U, N> const &b) {
  using std::declval;
  using R = decltype(declval<T>() * declval<U>());

  return std::inner_product(std::begin(a), std::end(a), std::begin(b), R{});
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

/**
 * @brief Meta function to turns a Vector<1, T> into T.
 */
template <typename T> struct decay_to_scalar {};
template <typename T, size_t N> struct decay_to_scalar<Vector<T, N>> {
  using type = Vector<T, N>;
};

template <typename T> struct decay_to_scalar<Vector<T, 1>> { using type = T; };
} // namespace Utils

#endif
