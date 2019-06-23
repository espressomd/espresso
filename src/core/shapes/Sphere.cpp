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

#include "Sphere.hpp"

#include <utils/math/sqr.hpp>

#include <cmath>

using namespace std;

namespace Shapes {
void Sphere::calculate_dist(const Utils::Vector3d &pos, double *dist,
                            double *vec) const {
  int i;
  double fac, c_dist;

  c_dist = 0.0;
  for (i = 0; i < 3; i++) {
    vec[i] = m_pos[i] - pos[i];
    c_dist += Utils::sqr(vec[i]);
  }
  c_dist = sqrt(c_dist);

  if (m_direction == -1) {
    /* apply force towards inside the sphere */
    *dist = m_rad - c_dist;
    fac = *dist / c_dist;
    for (i = 0; i < 3; i++)
      vec[i] *= fac;
  } else {
    /* apply force towards outside the sphere */
    *dist = c_dist - m_rad;
    fac = *dist / c_dist;
    for (i = 0; i < 3; i++)
      vec[i] *= -fac;
  }
}
} // namespace Shapes
