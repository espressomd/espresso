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
#include <initializer_list>
#include <iterator>
#include <vector>

#include <boost/serialization/array.hpp>

template <size_t n, typename Scalar> class Vector {
private:
  std::array<Scalar, n> d;

public:
  typedef typename std::array<Scalar, n>::iterator iterator;
  typedef typename std::array<Scalar, n>::const_iterator const_iterator;
  typedef typename std::array<Scalar, n>::reference reference;
  typedef typename std::array<Scalar, n>::const_reference const_reference;

  Vector() {}

  template <typename Container> explicit Vector(Container const &v) {
    assert(std::distance(std::begin(v), std::end(v)) <= n);
    std::copy(std::begin(v), std::end(v), d.begin());
  }

  template <typename InputIterator>
  Vector(InputIterator const &begin, InputIterator const &end) {
    assert(std::distance(begin, end) <= n);
    std::copy(begin, end, d.begin());
  }

  Vector(std::initializer_list<Scalar> l) {
    assert(l.size() <= n);
    std::copy(std::begin(l), std::end(l), d.begin());
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

  operator std::array<Scalar, n> const &() const { return d; }
  operator std::vector<Scalar>() const {
    return std::vector<Scalar>(std::begin(d), std::end(d));
  }

  std::vector<Scalar> as_vector() const {
    return std::vector<Scalar>(std::begin(d), std::end(d));
  }

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
    return cross(this, a);
  }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &d;
  }
};

typedef Vector<3, double> Vector3d;
typedef Vector<2, double> Vector2d;

#endif
