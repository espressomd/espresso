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

#include "Ellipsoid.hpp"
#include <cmath>
#include <iostream>

using namespace std;

#define SQR(A) ((A) * (A))

namespace Shapes {
int Ellipsoid::calculate_dist(const double *ppos, double *dist,
                              double *vec) const {

  Vector3d ppos_e;
  double l0;
  double l = 0.;
  double eps = 10.;
  int step = 0;
  double distance = 0.;

  // get particle position in reference frame of ellipsoid
  for (int i = 0; i < 3; i++) {
    ppos_e[i] = ppos[i] - m_pos[i];
  }

  // set appropriate initial point for Newton's method
  if (not inside_ellipsoid(ppos_e)) {
    l = *std::max_element(m_semiaxes.begin(), m_semiaxes.end()) * ppos_e.norm();
  }

  // find root via Newton's method
  while ((eps >= 1e-12) and (step < 10000)) {
    l0 = l;
    l -= lagrangian(ppos_e, l0) / lagrangian_derivative(ppos_e, l0);
    eps = std::abs(l - l0);
    step++;
  }

  for (int i = 0; i < 3; i++) {
    vec[i] =
        -m_direction *
        (SQR(m_semiaxes[i]) * ppos_e[i] / (l + SQR(m_semiaxes[i])) - ppos_e[i]);
    distance += SQR(vec[i]);
  }

  *dist = sqrt(distance);

  return 0;
}

bool Ellipsoid::inside_ellipsoid(const Vector3d ppos) const {
  bool is_inside = false;
  if (SQR(ppos[0] / m_semiaxes[0]) + SQR(ppos[1] / m_semiaxes[1]) +
          SQR(ppos[2] / m_semiaxes[2]) <=
      1)
    is_inside = true;
  return is_inside;
}

double Ellipsoid::lagrangian(const Vector3d ppos, const double l) const {
  return SQR(m_semiaxes[0]) * SQR(ppos[0]) * SQR(l + SQR(m_semiaxes[1])) *
             SQR(l + SQR(m_semiaxes[2])) +
         SQR(m_semiaxes[1]) * SQR(ppos[1]) * SQR(l + SQR(m_semiaxes[0])) *
             SQR(l + SQR(m_semiaxes[2])) +
         SQR(m_semiaxes[2]) * SQR(ppos[2]) * SQR(l + SQR(m_semiaxes[0])) *
             SQR(l + SQR(m_semiaxes[1])) -
         SQR(l + SQR(m_semiaxes[0])) * SQR(l + SQR(m_semiaxes[1])) *
             SQR(l + SQR(m_semiaxes[2]));
}

double Ellipsoid::lagrangian_derivative(const Vector3d ppos,
                                        const double l) const {
  return 2.0 * (SQR(m_semiaxes[0]) * SQR(ppos[0]) * (l + SQR(m_semiaxes[1])) *
                    SQR(l + SQR(m_semiaxes[2])) +
                SQR(m_semiaxes[0]) * SQR(ppos[0]) *
                    SQR(l + SQR(m_semiaxes[1])) * (l + SQR(m_semiaxes[2])) +
                SQR(m_semiaxes[1]) * SQR(ppos[1]) * (l + SQR(m_semiaxes[0])) *
                    SQR(l + SQR(m_semiaxes[2])) +
                SQR(m_semiaxes[1]) * SQR(ppos[1]) *
                    SQR(l + SQR(m_semiaxes[0])) * (l + SQR(m_semiaxes[2])) +
                SQR(m_semiaxes[2]) * SQR(ppos[2]) * (l + SQR(m_semiaxes[0])) *
                    SQR(l + SQR(m_semiaxes[1])) +
                SQR(m_semiaxes[2]) * SQR(ppos[2]) *
                    SQR(l + SQR(m_semiaxes[0])) * (l + SQR(m_semiaxes[1])) -
                (l + SQR(m_semiaxes[0])) * SQR(l + SQR(m_semiaxes[1])) *
                    SQR(l + SQR(m_semiaxes[2])) -
                SQR(l + SQR(m_semiaxes[0])) * (l + SQR(m_semiaxes[1])) *
                    SQR(l + SQR(m_semiaxes[2])) -
                SQR(l + SQR(m_semiaxes[0])) * SQR(l + SQR(m_semiaxes[1])) *
                    (l + SQR(m_semiaxes[2])));
}

} // namespace Shapes
