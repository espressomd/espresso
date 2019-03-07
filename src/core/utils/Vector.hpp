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

template <typename Scalar, std::size_t n> class Vector : public Utils::Array<Scalar, n> {
public:
  using Utils::Array<Scalar, n>::at;
  using Utils::Array<Scalar, n>::operator[];
  using Utils::Array<Scalar, n>::front;
  using Utils::Array<Scalar, n>::back;
  using Utils::Array<Scalar, n>::data;
  using Utils::Array<Scalar, n>::begin;
  using Utils::Array<Scalar, n>::cbegin;
  using Utils::Array<Scalar, n>::end;
  using Utils::Array<Scalar, n>::cend;
  using Utils::Array<Scalar, n>::empty;
  using Utils::Array<Scalar, n>::size;
  using Utils::Array<Scalar, n>::max_size;
  using Utils::Array<Scalar, n>::fill;
  using Utils::Array<Scalar, n>::broadcast;
  Vector() = default;
  Vector(Vector const &) = default;
  Vector &operator=(Vector const &) = default;

  void swap(Vector &rhs) { std::swap_ranges(begin(), end(), rhs.begin()); }

  template <typename Container>
  explicit Vector(Container const &v) : Vector(std::begin(v), std::end(v)) {}

  explicit Vector(Scalar const (&v)[n]) {
    std::copy_n(std::begin(v), n, begin());
  }

  Vector(std::initializer_list<Scalar> v)
      : Vector(std::begin(v), std::end(v)) {}

  template <typename InputIterator>
  Vector(InputIterator first, InputIterator last) {
    if (std::distance(first, last) == n) {
      std::copy_n(first, n, begin());
    } else {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }
  }

  /**
   * @brief Create a vector that has all entries set to
   *         one value.
   */
  static Vector<Scalar, n> broadcast(const Scalar &s) {
    Vector<Scalar, n> ret;
    std::fill(ret.begin(), ret.end(), s);

    return ret;
  }

  std::vector<Scalar> as_vector() const {
    return std::vector<Scalar>(begin(), end());
  }

  operator std::vector<Scalar>() const { return as_vector(); }

  inline Scalar dot(const Vector<Scalar, n> &b) const { return *this * b; }

  inline Scalar norm2() const { return (*this) * (*this); }
  inline Scalar norm() const { return sqrt(norm2()); }

  inline Vector &normalize(void) {
    const auto N = norm();
    if (N > Scalar(0)) {
      for (int i = 0; i < n; i++)
        this->operator[](i) /= N;
    }

    return *this;
  }

  static void cross(const Vector<Scalar, 3> &a, const Vector<Scalar, 3> &b,
                    Vector<Scalar, 3> &c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }

  static Vector<Scalar, 3> cross(const Vector<Scalar, 3> &a,
                                 const Vector<Scalar, 3> &b) {
    Vector<Scalar, 3> c;
    cross(a, b, c);
    return c;
  }

  inline Vector<Scalar, 3> cross(const Vector<Scalar, 3> &a) const {
    return cross(*this, a);
  }
};

// Useful typedefs

template <size_t N> using VectorXd = Vector<double, N>;
using Vector2d = VectorXd<2>;
using Vector3d = VectorXd<3>;
using Vector4d = VectorXd<4>;
using Vector6d = VectorXd<6>;
using Vector19d = VectorXd<19>;

using Vector3i = Vector<int, 3>;

namespace detail {
template <size_t N, typename T, typename Op>
Vector<T, N> binary_op(Vector<T, N> const &a, Vector<T, N> const &b, Op op) {
  Vector<T, N> ret;

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
bool all_of(Vector<T, N> const &a, Vector<T, N> const &b, Op op) {
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
bool operator<(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::less<T>());
}

template <size_t N, typename T>
bool operator>(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::greater<T>());
}

template <size_t N, typename T>
bool operator<=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::less_equal<T>());
}

template <size_t N, typename T>
bool operator>=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::greater_equal<T>());
}

template <size_t N, typename T>
bool operator==(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::all_of(a, b, std::equal_to<T>());
}

template <size_t N, typename T>
bool operator!=(Vector<T, N> const &a, Vector<T, N> const &b) {
  return not(a == b);
}

template <size_t N, typename T>
Vector<T, N> operator+(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::binary_op(a, b, std::plus<T>());
}

template <size_t N, typename T>
Vector<T, N> &operator+=(Vector<T, N> &a, Vector<T, N> const &b) {
  return detail::binary_op_assign(a, b, std::plus<T>());
}

template <size_t N, typename T>
Vector<T, N> operator-(Vector<T, N> const &a, Vector<T, N> const &b) {
  return detail::binary_op(a, b, std::minus<T>());
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
template <size_t N, typename T>
Vector<T, N> operator*(T const &a, Vector<T, N> const &b) {
  Vector<T, N> ret;

  std::transform(std::begin(b), std::end(b), std::begin(ret),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T>
Vector<T, N> operator*(Vector<T, N> const &b, T const &a) {
  Vector<T, N> ret;

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
template <size_t N, typename T>
T operator*(Vector<T, N> const &a, Vector<T, N> const &b) {
  return std::inner_product(std::begin(a), std::end(a), std::begin(b), T{});
}

/* Componentwise square root */
template <size_t N, typename T> Vector<T, N> sqrt(Vector<T, N> const &a) {
  using std::sqrt;
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), ret.begin(),
                 [](T const &v) { return sqrt(v); });

  return ret;
}

/**
 * @brief Meta function to turns a Vector<1, T> into T.
 */
template <typename T> struct decay_to_scalar {};
template <typename T, size_t N> struct decay_to_scalar<Vector<T, N>> {
  using type = Vector<T, N>;
};

template <typename T> struct decay_to_scalar<Vector<T, 1>> { using type = T; };

#endif
