/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
#ifndef UTILS_HPP
#define UTILS_HPP
/** \file utils.hpp
 *    Small functions that are useful not only for one modul.

 *  just some nice utilities...
 *  and some constants...
 *
*/

#include "Vector.hpp"
#include "utils/constants.hpp"
#include "utils/math/sqr.hpp"

#include <cassert>
#include <cmath>

/*************************************************************/
/** \name Vector and matrix operations for three dimensons.  */
/*************************************************************/
/*@{*/

/** calculates the scalar product of two vectors a nd b */
template <typename T1, typename T2> double scalar(const T1 &a, const T2 &b) {
  double d2 = 0.0;
  for (int i = 0; i < 3; i++)
    d2 += a[i] * b[i];
  return d2;
}

/** calculates the squared length of a vector */
template <typename T> double sqrlen(T const &v) { return scalar(v, v); }

/** calculates the vector product c of two vectors a and b */
template <typename T, typename U, typename V>
inline void vector_product(T const &a, U const &b, V &c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return;
}

/*@}*/

/*************************************************************/
/** \name Three dimensional grid operations                  */
/*************************************************************/
/*@{*/

/** get the linear index from the position (a,b,c) in a 3D grid
 *  of dimensions adim[]. returns linear index.
 *
 * @return        the linear index
 * @param a       x position
 * @param b       y position
 * @param c       z position
 * @param adim    dimensions of the underlying grid
 */
inline int get_linear_index(int a, int b, int c, const Vector3i &adim) {
  assert((a >= 0) && (a < adim[0]));
  assert((b >= 0) && (b < adim[1]));
  assert((c >= 0) && (c < adim[2]));

  return (a + adim[0] * (b + adim[1] * c));
}

/** get the position a[] from the linear index in a 3D grid
 *  of dimensions adim[].
 *
 * @param i       linear index
 * @param a       x position (return value)
 * @param b       y position (return value)
 * @param c       z position (return value)
 * @param adim    dimensions of the underlying grid
 */
inline void get_grid_pos(int i, int *a, int *b, int *c, const Vector3i &adim) {
  *a = i % adim[0];
  i /= adim[0];
  *b = i % adim[1];
  i /= adim[1];
  *c = i;
}

/*@}*/

/*************************************************************/
/** \name Distance calculations.  */
/*************************************************************/
/*@{*/

/** returns the distance between two position.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 */
inline double distance2(const Vector3d &a, const Vector3d &b) {
  return (a - b).norm2();
}

/*@}*/

#endif
