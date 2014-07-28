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
#include <algorithm>

template<int n, typename Scalar>
class Vector {
private:
  Scalar *d;
public:
  Vector() : d(new Scalar[n]) {};
  Vector(Scalar *a) : d(new Scalar[n]) {
    for(int i = 0; i < n; i++)
      d[i] = a[i];
  };
  ~Vector() { delete[] d; }
  Vector(const Vector& rhs)    : d(new Scalar[n]) {std::copy(rhs.d,rhs.d+n,d);}
  void swap(Vector& rhs)       {std::swap(d,rhs.d);}
  Vector& operator=(Vector& rhs) 
  {Vector tmp(rhs); swap(rhs); return *this;};
  Scalar &operator[](int i) {
    return d[i];
  };
  
};

typedef Vector<3, double> Vector3d;

