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
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <vector>

#include "utils/serialization/array.hpp"
#include <boost/serialization/access.hpp>

template <typename Scalar, std::size_t n> class Vector {
private:
  std::array<Scalar, n> d;

public:
  /* Concept Container requirements */
  using size_type = typename std::array<Scalar, n>::size_type;
  using difference_type = typename std::array<Scalar, n>::difference_type;
  using value_type = Scalar;
  using reference = Scalar &;
  using const_reference = const Scalar &;
  using iterator = typename std::array<Scalar, n>::iterator;
  using const_iterator = typename std::array<Scalar, n>::const_iterator;

  Vector() = default;
  Vector(Vector const &) = default;
  Vector &operator=(Vector const &) = default;

  void swap(Vector &rhs) { std::swap(d, rhs.d); }

  template <typename Container>
  explicit Vector(Container const &v) : Vector(std::begin(v), std::end(v)) {}

  explicit Vector(Scalar const (&v)[n]) {
    std::copy_n(std::begin(v), n, d.begin());
  }

  Vector(std::initializer_list<Scalar> v)
      : Vector(std::begin(v), std::end(v)) {}

  template <typename InputIterator>
  Vector(InputIterator begin, InputIterator end) {
    if (std::distance(begin, end) == n) {
      std::copy_n(begin, n, d.begin());
    } else {
      throw std::length_error(
          "Construction of Vector from Container of wrong length.");
    }
  }

  Scalar &operator[](int i) {
    assert(i < n);
    return d[i];
  }
  Scalar const &operator[](int i) const {
    assert(i < n);
    return d[i];
  }

  iterator begin() { return d.begin(); }
  const_iterator begin() const { return d.begin(); }
  const_iterator cbegin() const { return d.cbegin(); }

  iterator end() { return d.end(); }
  const_iterator end() const { return d.end(); }
  const_iterator cend() const { return d.cend(); }

  reference front() { return d.front(); }
  reference back() { return d.back(); }

  const_reference front() const { return d.front(); }
  const_reference back() const { return d.back(); }

  static constexpr size_t size() { return n; }
  Scalar const *data() const { return d.data(); }
  Scalar *data() { return d.data(); }
  size_type max_size() const { return d.max_size(); }
  bool empty() const { return d.empty(); }

  operator std::array<Scalar, n> const &() const { return d; }

  std::vector<Scalar> as_vector() const {
    return std::vector<Scalar>(std::begin(d), std::end(d));
  }

  operator std::vector<Scalar>() const { return as_vector(); }

  inline Scalar dot(const Vector<Scalar, n> &b) const { return *this * b; }

  inline Scalar norm2() const { return (*this) * (*this); }
  inline Scalar norm() const { return sqrt(norm2()); }

  inline Vector &normalize(void) {
    const auto N = norm();
    if (N > Scalar(0)) {
      for (int i = 0; i < n; i++)
        d[i] /= N;
    }

    return *this;
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

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &d;
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

  std::transform(a.begin(), a.end(), ret.begin(),
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

  std::transform(b.begin(), b.end(), ret.begin(),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T>
Vector<T, N> operator*(Vector<T, N> const &b, T const &a) {
  Vector<T, N> ret;

  std::transform(b.begin(), b.end(), ret.begin(),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T>
Vector<T, N> &operator*=(Vector<T, N> &b, T const &a) {
  std::transform(b.begin(), b.end(), b.begin(),
                 [a](T const &val) { return a * val; });
  return b;
}

/* Scalar division */
template <size_t N, typename T>
Vector<T, N> operator/(Vector<T, N> const &a, T const &b) {
  Vector<T, N> ret;

  std::transform(a.begin(), a.end(), ret.begin(),
                 [b](T const &val) { return val / b; });
  return ret;
}

template <size_t N, typename T>
Vector<T, N> &operator/=(Vector<T, N> &a, T const &b) {
  std::transform(a.begin(), a.end(), a.begin(),
                 [b](T const &val) { return val / b; });
  return a;
}

/* Scalar product */
template <size_t N, typename T>
T operator*(Vector<T, N> const &a, Vector<T, N> const &b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), T{});
}

/* Componentwise square root */
template <size_t N, typename T> Vector<T, N> sqrt(Vector<T, N> const &a) {
  using std::sqrt;
  Vector<T, N> ret;

  std::transform(a.begin(), a.end(), ret.begin(),
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
