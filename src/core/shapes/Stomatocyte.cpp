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

#include "Stomatocyte.hpp"
#include "utils.hpp"
#include <cmath>

using namespace std;

namespace Shapes {
void Stomatocyte::calculate_dist(const Utils::Vector3d &pos, double *dist,
                                 double *vec) const {

  using Utils::sqr;

  // Parameters

  int number;

  double mu, T0, T1, T1p, T2, T3, T3sqrt, T3p, T4, a, b, c, d, e, rad0, rad1,
      rad2, rad3, pt0x, pt0y, pt1x, pt1y, pt2x, pt2y, pt3x, pt3y, dst0, dst1,
      dst2, dst3, t0, t1, t2, t3, t4, ttota, distance, mindist, time0, time1,
      xd, yd, zd, xp, yp, zp, xpp, ypp, normal_x_3D, normal_y_3D, normal_z_3D;

  Utils::Vector3d closest_pos({-1.0, -1.0, -1.0});

  // Set the three dimensions of the stomatocyte

  a = m_outer_radius;
  b = m_inner_radius;
  c = m_layer_width;
  a = a * c;
  b = b * c;

  /***** Convert 3D coordinates to 2D planar coordinates *****/

  // Calculate the point on position + mu * orientation,
  // where the difference segment is orthogonal

  mu = (m_orientation * pos - m_position * m_orientation) /
       m_orientation.norm2();

  // Then the closest point to the line is

  closest_pos = m_position + mu * m_orientation;

  // So the shortest distance to the line is

  Utils::Vector2d dist_2D;
  dist_2D[0] = Utils::Vector3d(closest_pos - pos).norm();
  dist_2D[1] = mu * m_orientation.norm();

  /***** Use the obtained planar coordinates in distance function *****/

  // Calculate intermediate results which we need to determine
  // the distance, these are basically points where the different
  // segments of the 2D cut of the stomatocyte meet. The
  // geometry is rather complex however and details will not be given

  d = -sqrt(sqr(b) + 4 * b * c - 5 * sqr(c)) + a - 2 * c - b;

  e = (6 * sqr(c) - 2 * c * d + a * (-3 * c + d) +
       sqrt(-sqr(a - 2 * c) *
            (-2 * sqr(a) + 8 * a * c + sqr(c) + 6 * c * d + sqr(d)))) /
      (2 * (a - 2 * c));

  T0 = acos(
      (sqr(a) * d + 5 * sqr(c) * d + pow(d, 3) - sqr(a) * e - 5 * sqr(c) * e +
       6 * c * d * e - 3 * sqr(d) * e - 6 * c * sqr(e) + 4 * d * sqr(e) -
       2 * pow(e, 3) -
       (3 * c + e) *
           sqrt(-pow(a, 4) -
                sqr(5 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e + 2 * sqr(e)) +
                2 * sqr(a) *
                    (13 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e +
                     2 * sqr(e)))) /
      (2 * a * (9 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e + 2 * sqr(e))));

  T1p = acos(-(
      (39 * pow(c, 3) + 3 * c * sqr(d) + 31 * sqr(c) * e - 6 * c * d * e +
       sqr(d) * e + 12 * c * sqr(e) - 2 * d * sqr(e) + 2 * pow(e, 3) -
       sqr(a) * (3 * c + e) -
       d * sqrt(-pow(a, 4) -
                sqr(5 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e + 2 * sqr(e)) +
                2 * sqr(a) *
                    (13 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e +
                     2 * sqr(e))) +
       e * sqrt(-pow(a, 4) -
                sqr(5 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e + 2 * sqr(e)) +
                2 * sqr(a) *
                    (13 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e +
                     2 * sqr(e)))) /
      (4 * c * (9 * sqr(c) + sqr(d) + 6 * c * e - 2 * d * e + 2 * sqr(e)))));

  T1 = 3.0 * M_PI / 4.0 - T1p;

  T2 = e;

  T3sqrt = -sqr(b * c) *
           (sqr(a) + 9 * sqr(c) + 4 * c * d + d * (2 * b + d) -
            2 * a * (b + 2 * c + d)) *
           (sqr(a) + 9 * sqr(c) + 4 * c * d + sqr(d) - 2 * a * (b + 2 * c + d) +
            2 * b * (4 * c + d));

  T3sqrt = std::max(T3sqrt, 0.0);

  T3p = acos(
      -((-pow(a, 3) * b + 2 * pow(b, 3) * (2 * c + d) +
         3 * sqr(a) * b * (b + 2 * c + d) +
         sqr(b) * (25 * sqr(c) + 12 * c * d + 3 * sqr(d)) +
         b * (34 * pow(c, 3) + 25 * sqr(c) * d + 6 * c * sqr(d) + pow(d, 3)) -
         a * b *
             (2 * sqr(b) + 25 * sqr(c) + 12 * c * d + 3 * sqr(d) +
              6 * b * (2 * c + d)) +
         3 * sqrt(T3sqrt)) /
        (4 * b * c *
         (sqr(a) + sqr(b) + 13 * sqr(c) + 4 * c * d + sqr(d) +
          2 * b * (2 * c + d) - 2 * a * (b + 2 * c + d)))));

  T3 = 3.0 * M_PI / 4.0 - T3p;

  T4 = acos((b * (-a + b + 2 * c + d) *
                 (sqr(a) + 2 * sqr(b) + 9 * sqr(c) + 4 * c * d + sqr(d) +
                  2 * b * (2 * c + d) - 2 * a * (b + 2 * c + d)) -
             3 * sqrt(-sqr(b * c) *
                      (sqr(a) + 9 * sqr(c) + 4 * c * d + d * (2 * b + d) -
                       2 * a * (b + 2 * c + d)) *
                      (sqr(a) + 9 * sqr(c) + 4 * c * d + sqr(d) -
                       2 * a * (b + 2 * c + d) + 2 * b * (4 * c + d)))) /
            (2 * sqr(b) *
             (sqr(a) + sqr(b) + 13 * sqr(c) + 4 * c * d + sqr(d) +
              2 * b * (2 * c + d) - 2 * a * (b + 2 * c + d))));

  // Radii for the various parts of the swimmer

  rad0 = a;
  rad1 = 2.0 * c;
  rad2 = 2.0 * c;
  rad3 = b;

  // Center points for the circles

  Utils::Vector2d pt0, pt1, pt2, pt3;

  pt0 = {0.0, 0.0};
  pt1 = {3.0 * c + e, d - e};
  pt2 = {3.0 * c, d};
  pt3 = {0.0, a - b - 2 * c};

  // Distance of point of interest to center points

  dst0 = (pt0 - dist_2D).norm();
  dst1 = (pt1 - dist_2D).norm();
  dst2 = (pt2 - dist_2D).norm();
  dst3 = (pt3 - dist_2D).norm();

  // Now for the minimum distances, the fourth
  // is for the line segment

  double mdst[5];

  mdst[0] = (dst0 < rad0 ? rad0 - dst0 : dst0 - rad0);
  mdst[1] = (dst1 < rad1 ? rad1 - dst1 : dst1 - rad1);
  mdst[2] = (dst2 < rad2 ? rad2 - dst2 : dst2 - rad2);
  mdst[3] = (dst3 < rad3 ? rad3 - dst3 : dst3 - rad3);
  mdst[4] = fabs((-3.0 + 2.0 * sqrt(2.0)) * c - d + dist_2D[0] + dist_2D[1]) /
            sqrt(2.0);

  // Now determine at which time during the parametrization
  // the minimum distance to each segment is achieved

  ttota = T0 + T1 + T2 + T3 + T4;

  t0 = acos((dist_2D[1] - pt0[1]) / (dist_2D - pt0).norm());
  t1 = acos((dist_2D[0] - pt1[0]) / (dist_2D - pt1).norm());
  t2 = acos((dist_2D[1] - pt2[1]) / (dist_2D - pt2).norm());
  t3 = acos((dist_2D[1] - pt3[1]) / (dist_2D - pt3).norm());
  t4 = (3.0 * c - d + 2.0 * e + 2.0 * T0 + 2.0 * T1 - dist_2D[0] + dist_2D[1]) /
       (2.0 * ttota);

  // Now we can use these times to check whether or not the
  // point where the shortest distance is found is on the
  // segment of the circles or line that contributes to the
  // stomatocyte contour

  time0 = (T0 + T1) / ttota;
  time1 = (T0 + T1 + T2) / ttota;

  int io[5];

  io[0] = (0.0 <= t0 && t0 <= T0 ? 1 : 0);
  io[1] =
      (T1p <= t1 && t1 <= 3.0 * M_PI / 4.0 && (dist_2D[1] <= d - e) ? 1 : 0);
  io[2] =
      (T3p <= t2 && t2 <= 3.0 * M_PI / 4.0 && dist_2D[0] <= 3.0 * c ? 1 : 0);
  io[3] = (0.0 <= t3 && t3 <= T4 ? 1 : 0);
  io[4] = (time0 <= t4 && t4 <= time1 ? 1 : 0);

  // Now we only need to consider those distances for which
  // the io# flag is set to 1

  number = -1;
  mindist = -1.0;

  for (int i = 0; i < 5; i++) {
    if ((io[i] == 1) and ((mindist < 0.0) or (mindist > mdst[i]))) {
      number = i;
      mindist = mdst[i];
    }
  }

  // Now we know the number corresponding to the boundary
  // to which the point is closest, we know the distance,
  // but we still need the normal

  distance = -mindist;
  Utils::Vector2d normal({-1.0, -1.0});

  switch (number) {
  case 0:
    distance = (dst0 < rad0 ? -mindist : mindist);
    normal = (dist_2D - pt0);
    break;
  case 1:
    distance = (dst1 < rad1 ? -mindist : mindist);
    normal = (dist_2D - pt1);
    break;
  case 2:
    distance = (dst2 < rad2 ? -mindist : mindist);
    normal = (dist_2D - pt2);
    break;
  case 3:
    distance = (dst3 < rad3 ? mindist : -mindist);
    normal = -(dist_2D - pt3);
    break;
  case 4:
    if ((a - b + c - 2.0 * sqrt(2.0) * c - sqrt((b - c) * (b + 5.0 * c)) -
         dist_2D[0] - dist_2D[1]) > 0)
      distance = mindist;
    break;
  }

  normal.normalize();

  /***** Convert 2D normal to 3D coordinates *****/

  // Now that we have the normal in 2D we need to make a final
  // transformation to get it in 3D. The minimum distance stays
  // the same though. We first get the normalized direction vector.

  xd = m_orientation[0] / m_orientation.norm();
  yd = m_orientation[1] / m_orientation.norm();
  zd = m_orientation[2] / m_orientation.norm();

  // We now establish the rotation matrix required to go
  // form {0,0,1} to {xd,yd,zd}

  double matrix[9];

  if (xd * xd + yd * yd > 1.e-10) {
    // Rotation matrix

    matrix[0] = (yd * yd + xd * xd * zd) / (xd * xd + yd * yd);
    matrix[1] = (xd * yd * (zd - 1.0)) / (xd * xd + yd * yd);
    matrix[2] = xd;

    matrix[3] = (xd * yd * (zd - 1.0)) / (xd * xd + yd * yd);
    matrix[4] = (xd * xd + yd * yd * zd) / (xd * xd + yd * yd);
    matrix[5] = yd;

    matrix[6] = -xd;
    matrix[7] = -yd;
    matrix[8] = zd;
  } else {
    // The matrix is the identity matrix or reverses
    // or does a 180 degree rotation to take z -> -z

    matrix[0] = 1.0;
    matrix[1] = 0.0;
    matrix[2] = 0.0;

    matrix[3] = 0.0;
    matrix[4] = 1.0;
    matrix[5] = 0.0;

    matrix[6] = 0.0;
    matrix[7] = 0.0;

    if (zd > 0)
      matrix[8] = 1.0;
    else
      matrix[8] = -1.0;
  }

  // Next we determine the 3D vector between the center
  // of the stomatocyte and the point of interest

  Utils::Vector3d p(pos - m_position);

  // Now we use the inverse matrix to find the
  // position of the point with respect to the origin
  // of the z-axis oriented stomatocyte located
  // in the origin

  Utils::Vector2d pp;
  pp[0] = matrix[0] * p[0] + matrix[3] * p[1] + matrix[6] * p[2];
  pp[1] = matrix[1] * p[0] + matrix[4] * p[1] + matrix[7] * p[2];

  // Now use this direction to orient the normal

  if (pp.norm2() > 1.e-10) {
    // The point is off the rotational symmetry
    // axis of the stomatocyte

    pp.normalize();

    normal_x_3D = pp[0] * normal[0];
    normal_y_3D = pp[1] * normal[0];
    normal_z_3D = normal[1];
  } else {
    // The point is on the rotational symmetry
    // axis of the stomatocyte; a finite distance
    // away from the center the normal might have
    // an x and y component!

    normal_x_3D = 0.0;
    normal_y_3D = 0.0;
    normal_z_3D = (normal[1] > 0.0 ? 1.0 : -1.0);
  }

  // Now we need to transform the normal back to
  // the real coordinate system

  Utils::Vector3d normal_3D;

  normal_3D[0] = matrix[0] * normal_x_3D + matrix[1] * normal_y_3D +
                 matrix[2] * normal_z_3D;
  normal_3D[1] = matrix[3] * normal_x_3D + matrix[4] * normal_y_3D +
                 matrix[5] * normal_z_3D;
  normal_3D[2] = matrix[6] * normal_x_3D + matrix[7] * normal_y_3D +
                 matrix[8] * normal_z_3D;

  // Pass the values we obtained to ESPResSo

  for (int i = 0; i < 3; i++) {
    vec[i] = normal_3D[i] * distance;
  }

  *dist = std::copysign(distance, m_direction);

  // And we are done with the stomatocyte
}
} // namespace Shapes
