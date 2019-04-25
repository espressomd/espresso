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
/** \file
 *  Convenience functions for common operations on vectors.
 */

#include "utils/Vector.hpp"
#include "utils/constants.hpp"
#include "utils/math/sqr.hpp"

#include <boost/range/numeric.hpp>

/**************************************************************/
/** \name Vector and matrix operations for three dimensions.  */
/**************************************************************/
/*@{*/

/** calculates the scalar product of two vectors @p a and @p b */
template <typename T1, typename T2> double scalar(const T1 &a, const T2 &b) {
  double d2 = 0.0;
  for (int i = 0; i < 3; i++)
    d2 += a[i] * b[i];
  return d2;
}

/** calculates the squared length of a vector */
template <class T, size_t N> T sqrlen(const T (&v)[N]) {
    return boost::inner_product(v, v, T{});
}

/*@}*/

#endif
