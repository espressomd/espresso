/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#include "SpheroCylinder.hpp"
#include "utils.hpp"
#include <cmath>

using namespace std;

namespace Shapes {
int SpheroCylinder::calculate_dist(const double *ppos, double *dist,
                                   double *vec) const {
  int i;
  double d = 0.0;
  double ppos_local[3];

  for (i = 0; i < 3; i++) {
    ppos_local[i] = ppos[i] - m_pos[i];
    d += ppos_local[i] * m_axis[i];
  }

  if (abs(d) >= 0.5*m_length) {
    *dist = 0.0;

    for (i = 0; i < 3; i++) {
      vec[i] = ppos_local[i] - 0.5*m_length * m_axis[i] * Utils::sgn<double>(d);
      *dist += vec[i] * vec[i];
    }

    *dist = sqrt(*dist);

    if (*dist != 0.0)
      for (i = 0; i < 3; i++)
        vec[i] /= *dist;

    *dist -= m_rad;

    for (i = 0; i < 3; i++)
      vec[i] *= *dist;

    *dist *= m_direction;
  } else
    Cylinder::calculate_dist(ppos, dist, vec);
  return 0;
}
}
