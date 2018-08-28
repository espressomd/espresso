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

#include "Maze.hpp"
/* For box_l */
#include "grid.hpp"

using namespace std;

namespace Shapes {
int Maze::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int i, min_axis, cursph[3];
  double diasph, fac, c_dist, sph_dist, cyl_dist, temp_dis;
  double sph_vec[3], cyl_vec[3];

  diasph = box_l[0] / m_nsphere;

  /* First determine the distance to the sphere */
  c_dist = 0.0;
  for (i = 0; i < 3; i++) {
    cursph[i] = (int)(ppos[i] / diasph);
    sph_vec[i] = (cursph[i] + 0.5) * diasph - (ppos[i]);
    c_dist += Utils::sqr(sph_vec[i]);
  }
  c_dist = sqrt(c_dist);
  sph_dist = m_sphrad - c_dist;
  fac = sph_dist / c_dist;
  for (i = 0; i < 3; i++)
    cyl_vec[i] = sph_vec[i];
  for (i = 0; i < 3; i++)
    sph_vec[i] *= fac;

  /* Now calculate the cylinder stuff */
  /* have to check all 3 cylinders */
  min_axis = 2;
  cyl_dist = sqrt(cyl_vec[0] * cyl_vec[0] + cyl_vec[1] * cyl_vec[1]);

  if (m_dim > 0) {
    temp_dis = sqrt(cyl_vec[0] * cyl_vec[0] + cyl_vec[2] * cyl_vec[2]);
    if (temp_dis < cyl_dist) {
      min_axis = 1;
      cyl_dist = temp_dis;
    }

    if (m_dim > 1) {
      temp_dis = sqrt(cyl_vec[1] * cyl_vec[1] + cyl_vec[2] * cyl_vec[2]);
      if (temp_dis < cyl_dist) {
        min_axis = 0;
        cyl_dist = temp_dis;
      }
    }
  }
  cyl_vec[min_axis] = 0.;

  c_dist = cyl_dist;
  cyl_dist = m_cylrad - c_dist;
  fac = cyl_dist / c_dist;
  for (i = 0; i < 3; i++)
    cyl_vec[i] *= fac;

  /* Now decide between cylinder and sphere */
  if (sph_dist > 0) {
    if (sph_dist > cyl_dist) {
      *dist = sph_dist;
      for (i = 0; i < 3; i++)
        vec[i] = sph_vec[i];
    } else {
      *dist = cyl_dist;
      for (i = 0; i < 3; i++)
        vec[i] = cyl_vec[i];
    }
  } else {
    *dist = cyl_dist;
    for (i = 0; i < 3; i++)
      vec[i] = cyl_vec[i];
  }
  return 0;
}
} // namespace Shapes
