/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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
#include <vector>
#include <numeric>

#include "utils/serialization/array.hpp"

template <size_t n, typename Scalar> class Vector {
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

  template <typename Container> explicit Vector(Container v) {
    assert(std::distance(std::begin(v), std::end(v)) == n);
    std::copy(std::begin(v), std::end(v), d.begin());
  }

  template <typename T> explicit Vector(T const(&v)[n]) {
    std::copy(std::begin(v), std::end(v), d.begin());
  }

  Vector(std::initializer_list<Scalar> v) {
    /* Convert to static_assert in C++14 */
    assert(v.size() == n);
    std::copy(std::begin(v), std::end(v), d.begin());
  }

  template <typename InputIterator>
  Vector(InputIterator begin, InputIterator end) {
    assert(std::distance(begin, end) == n);
    std::copy(begin, end, d.begin());
  }

  Scalar &operator[](int i) { return d[i]; }
  Scalar const &operator[](int i) const { return d[i]; }

  iterator begin() { return d.begin(); }
  const_iterator begin() const { return d.begin(); }

  iterator end() { return d.end(); }
  const_iterator end() const { return d.end(); }

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

  inline Scalar dot(const Vector<n, Scalar> &b) const {
    Scalar sum = 0;
    for (int i = 0; i < n; i++)
      sum += d[i] * b[i];
    return sum;
  }

  inline Scalar norm2(void) const { return dot(*this); }

  inline Scalar norm(void) const { return sqrt(norm2()); }

  inline void normalize(void) {
    const auto N = norm();
    if (N > Scalar(0)) {
      for (int i = 0; i < n; i++)
        d[i] /= N;
    }
  }

  static void cross(const Vector<3, Scalar> &a, const Vector<3, Scalar> &b,
                    Vector<3, Scalar> &c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }

  static Vector<3, Scalar> cross(const Vector<3, Scalar> &a,
                                 const Vector<3, Scalar> &b) {
    Vector<3, Scalar> c;
    cross(a, b, c);
    return c;
  }

  inline Vector<3, Scalar> cross(const Vector<3, Scalar> &a) const {
    return cross(*this, a);
  }

  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &d;
  }
};

// Useful typedefs

typedef Vector<3, double> Vector3d;
typedef Vector<2, double> Vector2d;

namespace detail {
template <size_t N, typename T, typename Op>
Vector<N, T> binary_op(Vector<N, T> const &a, Vector<N, T> const &b, Op op) {
  Vector<N, T> ret;

  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(ret),
                 op);

  return ret;
}

template <size_t N, typename T, typename Op>
void binary_op_assign(Vector<N, T> &a, Vector<N, T> const &b, Op op) {
  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(a), op);
}

template <size_t N, typename T, typename Op>
bool all_of(Vector<N, T> const &a, Vector<N, T> const &b, Op op) {
  for (int i = 0; i < a.size(); i++) {
    /* Short circuit */
    if (!static_cast<bool>(op(a[i], b[i]))) {
      return false;
    }
  }

  return true;
}
}

template <size_t N, typename T>
bool operator<(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::all_of(a, b, std::less<T>());
}

template <size_t N, typename T>
bool operator>(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::all_of(a, b, std::greater<T>());
}

template <size_t N, typename T>
bool operator<=(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::all_of(a, b, std::less_equal<T>());
}

template <size_t N, typename T>
bool operator>=(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::all_of(a, b, std::greater_equal<T>());
}

template <size_t N, typename T>
bool operator==(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::all_of(a, b, std::equal_to<T>());
}

template <size_t N, typename T>
bool operator!=(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::all_of(a, b, std::not_equal_to<T>());
}

template <size_t N, typename T>
Vector<N, T> operator+(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::binary_op(a, b, std::plus<T>());
}

template <size_t N, typename T>
void operator+=(Vector<N, T> &a, Vector<N, T> const &b) {
  return detail::binary_op_assign(a, b, std::plus<T>());
}

template <size_t N, typename T>
Vector<N, T> operator-(Vector<N, T> const &a, Vector<N, T> const &b) {
  return detail::binary_op(a, b, std::minus<T>());
}

template <size_t N, typename T> Vector<N, T> operator-(Vector<N, T> const &a) {
  Vector<N, T> ret;

  std::transform(a.begin(), a.end(), ret.begin(), [](T const &v) { return -v; });

  return ret;
}

template <size_t N, typename T>
void operator-=(Vector<N, T> &a, Vector<N, T> const &b) {
  return detail::binary_op_assign(a, b, std::minus<T>());
}

/* Scalar multiplication */
template <size_t N, typename T>
Vector<N, T> operator*(T const &a, Vector<N, T> const &b) {
  Vector<N, T> ret;

  std::transform(b.begin(), b.end(), ret.begin(),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T>
Vector<N, T> operator*(Vector<N, T> const &b, T const &a) {
  Vector<N, T> ret;

  std::transform(b.begin(), b.end(), ret.begin(),
                 [a](T const &val) { return a * val; });

  return ret;
}

template <size_t N, typename T> void operator*=(Vector<N, T> &b, T const &a) {
  std::transform(b.begin(), b.end(), b.begin(),
                 [a](T const &val) { return a * val; });
}

/* Scalar division */
template <size_t N, typename T>
Vector<N, T> operator/(Vector<N, T> const &a, T const &b) {
  Vector<N, T> ret;

  std::transform(a.begin(), a.end(), ret.begin(),
                 [b](T const &val) { return val / b; });
  return ret;
}

template <size_t N, typename T> void operator/=(Vector<N, T> &a, T const &b) {
  std::transform(a.begin(), a.end(), a.begin(),
                 [b](T const &val) { return val / b; });
}

/* Scalar product */
template <size_t N, typename T>
T operator*(Vector<N, T> const &a, Vector<N, T> const &b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), T{});
}

/* Componentwise square route */
template <size_t N, typename T> 
Vector<N, T> sqrt(Vector<N, T> const& a) {
  using std::sqrt;
  Vector<N, T> ret;

  std::transform(a.begin(), a.end(), ret.begin(), [](T const& v) {
      return sqrt(v);
    });

  return ret;
}

#endif
