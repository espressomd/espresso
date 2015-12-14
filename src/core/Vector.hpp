/*
  Copyright (C) 2014 The ESPResSo project
  
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

#ifndef __UTILS_VECTOR_HPP
#define __UTILS_VECTOR_HPP

#include <algorithm>
#include <cmath>
#include <vector>
#ifdef HAVE_CXX11
#include <initializer_list>
#endif

template<int n, typename Scalar>
class Vector {
 private:
  Scalar d[n];
  
 public:
  Vector() {}

  explicit Vector(const Scalar * const a) {
    for (int i = 0; i < n; i++)
      d[i] = a[i];
  }

  Vector(const Vector& rhs) {
    std::copy(rhs.d, rhs.d+n, d);
  }

  Vector(const std::vector<Scalar> &rhs) {
    std::copy(rhs.d, rhs.d+n, d);
  }
  
#ifdef HAVE_CXX11
  Vector(std::initializer_list<Scalar> l) {
    std::copy(l.begin(), l.begin() + n, d);
  }
#endif
  
  Scalar *begin() const {
    return d;
  }

  Scalar *end() const {
    return d + n;
  }
 
  static int size() { return n; }
    
  Vector& operator=(Vector& rhs) {
    std::copy(rhs.d, rhs.d + n, d);
    return *this;
  }

  Vector& operator=(const std::vector<Scalar>& rhs) {
    std::copy(rhs.begin(), rhs.end(), d);
    return *this;
  }
  
  Scalar &operator[](int i) {
    return d[i];
  }

  Scalar const & operator[](int i) const {
    return d[i];
  }
  
  inline Scalar dot(const Vector<n, Scalar> &b) const {
    Scalar sum = 0;
    for (int i = 0; i < n; i++)
      sum += d[i]*b[i];
    return sum;
  }

  inline Scalar norm2(void) const {
    return dot(*this);
  }

  inline Scalar norm(void) const {
    return std::sqrt(norm2());
  }

  inline void normalize(void) {
    Scalar N = norm();
    if (norm() > 0) {
      for (int i = 0; i < n; i++)
        d[i] /= N;
    }
  }

  inline void cross(const Vector<3, Scalar> &a, const Vector<3, Scalar> &b, Vector<3, Scalar> &c) const {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
  }
  
  inline Vector<3, Scalar> cross(const Vector<3, Scalar> &a, const Vector<3, Scalar> &b) const {
    Vector<3, Scalar> c;
    cross(a,b,c);
    return c;
  }
  
  inline Vector<3, Scalar> cross(const Vector<3, Scalar> &a) const {
    return cross(this,a);
  }
};

typedef Vector<3, double> Vector3d;

#endif
