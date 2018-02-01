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

#include "Pore.hpp"
#include "utils.hpp"

#include <cmath>

using namespace std;


namespace Shapes {
int Pore::calculate_dist(const double *ppos, double *dist, double *vec) const {
  double c_dist[3]; /* cartesian distance from pore center */
  double z, r;      /* cylindrical coordinates, coordinate system parallel to
                       pore, origin at pore centera */
  double z_vec[3],
      r_vec[3]; /* cartesian vectors that correspond to these coordinates */
  double e_z[3], e_r[3]; /* unit vectors in the cylindrical coordinate system */
  /* helper variables, for performance reasons should the be move the the
   * constraint struct*/
  double z_left, z_right;
  /* and another helper that is hopefully optmized out */
  double norm;
  double c1_r, c1_z, c2_r, c2_z;
  double cone_vector_r, cone_vector_z, p1_r, p1_z, dist_vector_z, dist_vector_r,
      temp;

  auto const half_length = 0.5 * m_length;

  auto const slope =
      (m_rad_right - m_rad_left) / 2. / (half_length - m_smoothing_radius);
  auto const slope2 = (m_outer_rad_right - m_outer_rad_left) / 2. /
                      (half_length - m_smoothing_radius);

  /* compute the position relative to the center of the pore */
  for (int i = 0; i < 3; i++) {
    c_dist[i] = ppos[i] - m_pos[i];
  }

  /* compute the component parallel to the pore axis */
  z = 0.;
  for (int i = 0; i < 3; i++) {
    z += (c_dist[i] * m_axis[i]);
  }

  /* decompose the position into parallel and perpendicular to the axis */
  r = 0.;
  for (int i = 0; i < 3; i++) {
    z_vec[i] = z * m_axis[i];
    r_vec[i] = c_dist[i] - z_vec[i];
    r += r_vec[i] * r_vec[i];
  }
  r = sqrt(r);

  /* calculate norm and unit vectors for both */
  norm = 0;
  for (int i = 0; i < 3; i++)
    norm += z_vec[i] * z_vec[i];
  norm = sqrt(norm);
  for (int i = 0; i < 3; i++)
    e_z[i] = m_axis[i];
  norm = 0;
  for (int i = 0; i < 3; i++)
    norm += r_vec[i] * r_vec[i];
  norm = sqrt(norm);
  for (int i = 0; i < 3; i++)
    e_r[i] = r_vec[i] / norm;

  /* c?_r/z and are the centers of the circles that are used to smooth
   * the entrance of the pore in cylindrical coordinates*/
  c1_z = -(half_length - m_smoothing_radius);
  c2_z = +(half_length - m_smoothing_radius);
  z_left = c1_z -
           Utils::sgn<double>(slope) *
               sqrt(slope * slope / (1 + slope * slope)) * m_smoothing_radius;
  z_right = c2_z +
            Utils::sgn<double>(slope) *
                sqrt(slope * slope / (1 + slope * slope)) * m_smoothing_radius;

  c1_r = m_rad_left + slope * (z_left + half_length) +
         sqrt(m_smoothing_radius * m_smoothing_radius - Utils::sqr(z_left - c1_z));
  c2_r = m_rad_left + slope * (z_right + half_length) +
         sqrt(m_smoothing_radius * m_smoothing_radius - Utils::sqr(z_right - c2_z));
  c1_r = m_rad_left + m_smoothing_radius;
  c2_r = m_rad_right + m_smoothing_radius;

  double c1_or = m_outer_rad_left - m_smoothing_radius;
  double c2_or = m_outer_rad_right - m_smoothing_radius;

  /* Check if we are in the region of the left wall */
  if (((r >= c1_r) && (r <= c1_or) && (z <= c1_z))) {
    dist_vector_z = -z - half_length;
    dist_vector_r = 0;
    *dist = -z - half_length;
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }
  /* Check if we are in the region of the right wall */
  if (((r >= c2_r) && (r < c2_or) && (z >= c2_z))) {
    dist_vector_z = -z + half_length;
    dist_vector_r = 0;
    *dist = +z - half_length;
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }

  /* check if the particle should feel the smoothed ends or the middle of the
   * pore */
  /* calculate aufpunkt in z direction first.   */

  /* the distance of the particle from the pore cylinder/cone calculated by
   * projection on the
   * cone normal. Should be > 0 if particle is inside the pore */

  cone_vector_z = 1 / sqrt(1 + slope * slope);
  cone_vector_r = slope / sqrt(1 + slope * slope);

  double cone_vector_z_o = 1 / sqrt(1 + slope2 * slope2);
  double cone_vector_r_o = slope2 / sqrt(1 + slope2 * slope2);

  p1_r =
      c1_r +
      ((r - c1_r) * cone_vector_r + (z - c1_z) * cone_vector_z) * cone_vector_r;
  p1_z =
      c1_z +
      ((r - c1_r) * cone_vector_r + (z - c1_z) * cone_vector_z) * cone_vector_z;

  double p2_r = c1_or +
                ((r - c1_or) * cone_vector_r_o + (z - c1_z) * cone_vector_z_o) *
                    cone_vector_r_o;
  double p2_z = c1_z +
                ((r - c1_or) * cone_vector_r_o + (z - c1_z) * cone_vector_z_o) *
                    cone_vector_z_o;

  dist_vector_r = p1_r - r;
  dist_vector_z = p1_z - z;

  double dist_vector_r_o = p2_r - r;
  double dist_vector_z_o = p2_z - z;

  if (p1_z >= c1_z && p1_z <= c2_z && dist_vector_r >= 0) {
    temp = sqrt(dist_vector_r * dist_vector_r + dist_vector_z * dist_vector_z);
    *dist = temp - m_smoothing_radius;
    dist_vector_r -= dist_vector_r / temp * m_smoothing_radius;
    dist_vector_z -= dist_vector_z / temp * m_smoothing_radius;
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }

  if (p2_z >= c1_z && p2_z <= c2_z && dist_vector_r_o <= 0) {
    temp = sqrt(dist_vector_r_o * dist_vector_r_o +
                dist_vector_z_o * dist_vector_z_o);
    *dist = temp - m_smoothing_radius;
    dist_vector_r_o -= dist_vector_r_o / temp * m_smoothing_radius;
    dist_vector_z_o -= dist_vector_z_o / temp * m_smoothing_radius;
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r_o * e_r[i] - dist_vector_z_o * e_z[i];
    return 0;
  }

  /* Check if we are in the range of the left smoothing circle */
  if (p1_z <= c1_z && r <= c1_r) {
    /* distance from the smoothing center */
    norm = sqrt((z - c1_z) * (z - c1_z) + (r - c1_r) * (r - c1_r));
    *dist = norm - m_smoothing_radius;
    dist_vector_r = (m_smoothing_radius / norm - 1) * (r - c1_r);
    dist_vector_z = (m_smoothing_radius / norm - 1) * (z - c1_z);
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }
  /* upper left smoothing circle */
  if (p2_z <= c1_z && r >= c1_or) {
    /* distance from the smoothing center */
    norm = sqrt((z - c1_z) * (z - c1_z) + (r - c1_or) * (r - c1_or));
    *dist = norm - m_smoothing_radius;
    dist_vector_r = (m_smoothing_radius / norm - 1) * (r - c1_or);
    dist_vector_z = (m_smoothing_radius / norm - 1) * (z - c1_z);
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }
  /* Check if we are in the range of the right smoothing circle */
  if (p1_z >= c2_z && r <= c2_r) {
    norm = sqrt((z - c2_z) * (z - c2_z) + (r - c2_r) * (r - c2_r));
    *dist = norm - m_smoothing_radius;
    dist_vector_r = (m_smoothing_radius / norm - 1) * (r - c2_or);
    dist_vector_z = (m_smoothing_radius / norm - 1) * (z - c2_z);
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }
  /* Check if we are in the range of the upper right smoothing circle */
  if (p2_z >= c2_z && r >= c2_or) {
    norm = sqrt((z - c2_z) * (z - c2_z) + (r - c2_or) * (r - c2_or));
    *dist = norm - m_smoothing_radius;
    dist_vector_r = (m_smoothing_radius / norm - 1) * (r - c2_or);
    dist_vector_z = (m_smoothing_radius / norm - 1) * (z - c2_z);
    for (int i = 0; i < 3; i++)
      vec[i] = -dist_vector_r * e_r[i] - dist_vector_z * e_z[i];
    return 0;
  }
  *dist = -1e99;
  vec[0] = vec[1] = vec[2] = 1e99;
  return 0;
  //  exit(printf("should never be reached, z %f, r%f\n",z, r));
}
}
