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

#include "Cylinder.hpp"
#include <cmath>

using namespace std;

#define SQR(A) ((A) * (A))

namespace Shapes {
int Cylinder::calculate_dist(const double *ppos, double *dist,
                             double *vec) const {
  double d_per, d_par, d_real, d_per_vec[3], d_par_vec[3], d_real_vec[3];
  auto const half_length = 0.5 * m_length;

  d_real = 0.0;
  for (int i = 0; i < 3; i++) {
    d_real_vec[i] = ppos[i] - m_pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);

  d_par = 0.;
  for (int i = 0; i < 3; i++) {
    d_par += (d_real_vec[i] * m_axis[i]);
  }

  for (int i = 0; i < 3; i++) {
    d_par_vec[i] = d_par * m_axis[i];
    d_per_vec[i] = ppos[i] - (m_pos[i] + d_par_vec[i]);
  }

  d_per = sqrt(SQR(d_real) - SQR(d_par));
  d_par = fabs(d_par);

  if (m_direction == -1) {
    /* apply force towards inside cylinder */
    d_per = m_rad - d_per;
    d_par = half_length - d_par;
    if (d_per < d_par or m_open) {
      *dist = d_per;
      for (int i = 0; i < 3; i++) {
        vec[i] = -d_per_vec[i] * d_per / (m_rad - d_per);
      }
    } else {
      *dist = d_par;
      for (int i = 0; i < 3; i++) {
        vec[i] = -d_par_vec[i] * d_par / (half_length - d_par);
      }
    }
  } else {
    /* apply force towards outside cylinder */
    d_per = d_per - m_rad;
    d_par = d_par - half_length;
    if (d_par < 0) {
      *dist = d_per;
      for (int i = 0; i < 3; i++) {
        vec[i] = d_per_vec[i] * d_per / (d_per + m_rad);
      }
    } else if (d_per < 0) {
      *dist = d_par;
      for (int i = 0; i < 3; i++) {
        vec[i] = d_par_vec[i] * d_par / (d_par + half_length);
      }
    } else {
      *dist = sqrt(SQR(d_par) + SQR(d_per));
      for (int i = 0; i < 3; i++) {
        vec[i] = d_per_vec[i] * d_per / (d_per + m_rad) +
                 d_par_vec[i] * d_par / (d_par + half_length);
      }
    }
  }
  return 0;
}
}
