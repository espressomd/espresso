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

#include "Stomatocyte.hpp"
#include <cmath>

using namespace std;

namespace Shapes {
int Stomatocyte::calculate_dist(const double *ppos, double *dist,
                                double *vec) const {
  // Parameters

  int io0, io1, io2, io3, io4, number;

  double x_2D, y_2D, mu, T0, T1, T1p, T2, T3, T3sqrt, T3p, T4sqrt, T4, a, b, c,
      d, e, rad0, rad1, rad2, rad3, pt0x, pt0y, pt1x, pt1y, pt2x, pt2y, pt3x,
      pt3y, dst0, dst1, dst2, dst3, mdst0, mdst1, mdst2, mdst3, mdst4, t0, t1,
      t2, t3, t4, ttota, distance, mindist, normal_x, normal_y, time0, time1,
      xd, x, yd, y, zd, z, normal_3D_x, normal_3D_y, normal_3D_z, xp, yp, zp,
      xpp, ypp, normal_x_3D, normal_y_3D, normal_z_3D, sin_xy, cos_xy;

  double closest_point_3D[3] = {-1.0, -1.0, -1.0};

  // Set the three dimensions of the stomatocyte

  a = m_outer_radius;
  b = m_inner_radius;
  c = m_layer_width;
  a = a * c;
  b = b * c;

  // Set the position and orientation of the stomatocyte

  double stomatocyte_3D_position[3] = {m_position_x, m_position_y,
                                       m_position_z};

  double stomatocyte_3D_orientation[3] = {m_orientation_x, m_orientation_y,
                                          m_orientation_z};

  // Set the point for which we want to know the distance

  double point_3D[3];

  point_3D[0] = ppos[0];
  point_3D[1] = ppos[1];
  point_3D[2] = ppos[2];

  /***** Convert 3D coordinates to 2D planar coordinates *****/

  // Calculate the point on position + mu * orientation,
  // where the difference segment is orthogonal

  mu = (stomatocyte_3D_orientation[0] * point_3D[0] +
        stomatocyte_3D_orientation[1] * point_3D[1] +
        stomatocyte_3D_orientation[2] * point_3D[2] -
        stomatocyte_3D_position[0] * stomatocyte_3D_orientation[0] -
        stomatocyte_3D_position[1] * stomatocyte_3D_orientation[1] -
        stomatocyte_3D_position[2] * stomatocyte_3D_orientation[2]) /
       (stomatocyte_3D_orientation[0] * stomatocyte_3D_orientation[0] +
        stomatocyte_3D_orientation[1] * stomatocyte_3D_orientation[1] +
        stomatocyte_3D_orientation[2] * stomatocyte_3D_orientation[2]);

  // Then the closest point to the line is

  closest_point_3D[0] =
      stomatocyte_3D_position[0] + mu * stomatocyte_3D_orientation[0];
  closest_point_3D[1] =
      stomatocyte_3D_position[1] + mu * stomatocyte_3D_orientation[1];
  closest_point_3D[2] =
      stomatocyte_3D_position[2] + mu * stomatocyte_3D_orientation[2];

  // So the shortest distance to the line is

  x_2D = sqrt((closest_point_3D[0] - point_3D[0]) *
                  (closest_point_3D[0] - point_3D[0]) +
              (closest_point_3D[1] - point_3D[1]) *
                  (closest_point_3D[1] - point_3D[1]) +
              (closest_point_3D[2] - point_3D[2]) *
                  (closest_point_3D[2] - point_3D[2]));

  y_2D =
      mu * sqrt(stomatocyte_3D_orientation[0] * stomatocyte_3D_orientation[0] +
                stomatocyte_3D_orientation[1] * stomatocyte_3D_orientation[1] +
                stomatocyte_3D_orientation[2] * stomatocyte_3D_orientation[2]);

  /***** Use the obtained planar coordinates in distance function *****/

  // Calculate intermediate results which we need to determine
  // the distance, these are basically points where the different
  // segments of the 2D cut of the stomatocyte meet. The
  // geometry is rather complex however and details will not be given

  d = -sqrt(b * b + 4.0 * b * c - 5.0 * c * c) + a - 2.0 * c - b;

  e = (6.0 * c * c - 2.0 * c * d + a * (-3.0 * c + d) +
       sqrt(-(a - 2.0 * c) * (a - 2.0 * c) *
            (-2.0 * a * a + 8.0 * a * c + c * c + 6.0 * c * d + d * d))) /
      (2.0 * (a - 2.0 * c));

  T0 = acos(
      (a * a * d + 5.0 * c * c * d + d * d * d - a * a * e - 5.0 * c * c * e +
       6.0 * c * d * e - 3.0 * d * d * e - 6.0 * c * e * e + 4.0 * d * e * e -
       2.0 * e * e * e -
       (3.0 * c + e) * sqrt(-a * a * a * a -
                            (5.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e +
                             2.0 * e * e) *
                                (5.0 * c * c + d * d + 6.0 * c * e -
                                 2.0 * d * e + 2.0 * e * e) +
                            2.0 * a * a * (13.0 * c * c + d * d + 6.0 * c * e -
                                           2.0 * d * e + 2.0 * e * e))) /
      (2.0 * a *
       (9.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e + 2.0 * e * e)));

  T1p = acos(-(39.0 * c * c * c + 3.0 * c * d * d + 31.0 * c * c * e -
               6.0 * c * d * e + d * d * e + 12.0 * c * e * e -
               2.0 * d * e * e + 2.0 * e * e * e - a * a * (3.0 * c + e) -
               d * sqrt(-a * a * a * a -
                        (5.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e +
                         2.0 * e * e) *
                            (5.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e +
                             2.0 * e * e) +
                        2.0 * a * a * (13.0 * c * c + d * d + 6.0 * c * e -
                                       2.0 * d * e + 2.0 * e * e)) +
               e * sqrt(-a * a * a * a -
                        (5.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e +
                         2.0 * e * e) *
                            (5.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e +
                             2.0 * e * e) +
                        2.0 * a * a * (13.0 * c * c + d * d + 6.0 * c * e -
                                       2.0 * d * e + 2.0 * e * e))) /
             (4.0 * c *
              (9.0 * c * c + d * d + 6.0 * c * e - 2.0 * d * e + 2.0 * e * e)));

  T1 = 3.0 * M_PI / 4.0 - T1p;

  T2 = e;

  T3sqrt = -b * b * c * c *
           (a * a * a * a - 4.0 * a * a * a * (b + 2.0 * c + d) +
            4.0 * b * b * d * (4.0 * c + d) +
            (9.0 * c * c + 4.0 * c * d + d * d) *
                (9.0 * c * c + 4.0 * c * d + d * d) +
            4.0 * b * (18.0 * c * c * c + 17.0 * c * c * d + 6.0 * c * d * d +
                       d * d * d) +
            2.0 * a * a * (2.0 * b * b + 17.0 * c * c + 12.0 * c * d +
                           3.0 * d * d + 6.0 * b * (2.0 * c + d)) -
            4.0 * a * (18.0 * c * c * c + 17.0 * c * c * d + 6.0 * c * d * d +
                       d * d * d + 2.0 * b * b * (2.0 * c + d) +
                       b * (17.0 * c * c + 12.0 * c * d + 3.0 * d * d)));

  if (T3sqrt < 0.0)
    T3sqrt = 0.0;

  T3p = acos(
      -(-a * a * a * b + 4.0 * b * b * b * c + 25.0 * b * b * c * c +
        34.0 * b * c * c * c + 2.0 * b * b * b * d + 12.0 * b * b * c * d +
        25.0 * b * c * c * d + 3.0 * b * b * d * d + 6.0 * b * c * d * d +
        b * d * d * d + 3.0 * a * a * b * (b + 2.0 * c + d) -
        a * b * (2.0 * b * b + 25.0 * c * c + 12.0 * c * d + 3.0 * d * d +
                 6.0 * b * (2.0 * c + d)) +
        3.0 * sqrt(T3sqrt)) /
      (4.0 * b * c * (a * a + b * b + 13.0 * c * c + 4.0 * c * d + d * d +
                      2.0 * b * (2.0 * c + d) - 2.0 * a * (b + 2.0 * c + d))));

  T3 = 3.0 * M_PI / 4.0 - T3p;

  T4sqrt = -b * b * c * c *
           (a * a * a * a - 4.0 * a * a * a * (b + 2.0 * c + d) +
            4.0 * b * b * d * (4.0 * c + d) +
            (9.0 * c * c + 4.0 * c * d + d * d) *
                (9.0 * c * c + 4.0 * c * d + d * d) +
            4.0 * b * (18.0 * c * c * c + 17.0 * c * c * d + 6.0 * c * d * d +
                       d * d * d) +
            2.0 * a * a * (2.0 * b * b + 17.0 * c * c + 12.0 * c * d +
                           3.0 * d * d + 6.0 * b * (2.0 * c + d)) -
            4.0 * a * (18.0 * c * c * c + 17.0 * c * c * d + 6.0 * c * d * d +
                       d * d * d + 2.0 * b * b * (2.0 * c + d) +
                       b * (17.0 * c * c + 12.0 * c * d + 3.0 * d * d)));

  if (T4sqrt < 0.0)
    T4sqrt = 0.0;

  T4 = acos(
      (-a * a * a * b + 2.0 * b * b * b * b + 8.0 * b * b * b * c +
       17.0 * b * b * c * c + 18.0 * b * c * c * c + 4.0 * b * b * b * d +
       12.0 * b * b * c * d + 17.0 * b * c * c * d + 3.0 * b * b * d * d +
       6.0 * b * c * d * d + b * d * d * d +
       3.0 * a * a * b * (b + 2.0 * c + d) -
       a * b * (4.0 * b * b + 17.0 * c * c + 12.0 * c * d + 3.0 * d * d +
                6.0 * b * (2.0 * c + d)) -
       3.0 * sqrt(T4sqrt)) /
      (2.0 * b * b * (a * a + b * b + 13.0 * c * c + 4.0 * c * d + d * d +
                      2.0 * b * (2.0 * c + d) - 2.0 * a * (b + 2.0 * c + d))));

  // Radii for the various parts of the swimmer

  rad0 = a;
  rad1 = 2.0 * c;
  rad2 = 2.0 * c;
  rad3 = b;

  // Center points for the circles

  pt0x = 0.0;
  pt0y = 0.0;

  pt1x = 3.0 * c + e;
  pt1y = d - e;

  pt2x = 3.0 * c;
  pt2y = d;

  pt3x = 0.0;
  pt3y = a - b - 2 * c;

  // Distance of point of interest to center points

  dst0 = sqrt((pt0x - x_2D) * (pt0x - x_2D) + (pt0y - y_2D) * (pt0y - y_2D));
  dst1 = sqrt((pt1x - x_2D) * (pt1x - x_2D) + (pt1y - y_2D) * (pt1y - y_2D));
  dst2 = sqrt((pt2x - x_2D) * (pt2x - x_2D) + (pt2y - y_2D) * (pt2y - y_2D));
  dst3 = sqrt((pt3x - x_2D) * (pt3x - x_2D) + (pt3y - y_2D) * (pt3y - y_2D));

  // Now for the minimum distances, the fourth
  // is for the line segment

  mdst0 = (dst0 < rad0 ? rad0 - dst0 : dst0 - rad0);
  mdst1 = (dst1 < rad1 ? rad1 - dst1 : dst1 - rad1);
  mdst2 = (dst2 < rad2 ? rad2 - dst2 : dst2 - rad2);
  mdst3 = (dst3 < rad3 ? rad3 - dst3 : dst3 - rad3);

  mdst4 = fabs((-3.0 + 2.0 * sqrt(2.0)) * c - d + x_2D + y_2D) / sqrt(2.0);

  // Now determine at which time during the parametrization
  // the minimum distance to each segment is achieved

  ttota = T0 + T1 + T2 + T3 + T4;

  t0 = acos((y_2D - pt0y) / sqrt((x_2D - pt0x) * (x_2D - pt0x) +
                                 (y_2D - pt0y) * (y_2D - pt0y)));
  t1 = acos((x_2D - pt1x) / sqrt((x_2D - pt1x) * (x_2D - pt1x) +
                                 (y_2D - pt1y) * (y_2D - pt1y)));
  t2 = acos((y_2D - pt2y) / sqrt((x_2D - pt2x) * (x_2D - pt2x) +
                                 (y_2D - pt2y) * (y_2D - pt2y)));
  t3 = acos((y_2D - pt3y) / sqrt((x_2D - pt3x) * (x_2D - pt3x) +
                                 (y_2D - pt3y) * (y_2D - pt3y)));

  t4 = (3.0 * c - d + 2.0 * e + 2.0 * T0 + 2.0 * T1 - x_2D + y_2D) /
       (2.0 * ttota);

  // Now we can use these times to check whether or not the
  // point where the shortest distance is found is on the
  // segment of the circles or line that contributes to the
  // stomatocyte contour

  time0 = (T0 + T1) / ttota;
  time1 = (T0 + T1 + T2) / ttota;

  io0 = (0.0 <= t0 && t0 <= T0 ? 1 : 0);
  io1 = (T1p <= t1 && t1 <= 3.0 * M_PI / 4.0 && (y_2D <= d - e) ? 1 : 0);
  io2 = (T3p <= t2 && t2 <= 3.0 * M_PI / 4.0 && x_2D <= 3.0 * c ? 1 : 0);
  io3 = (0.0 <= t3 && t3 <= T4 ? 1 : 0);

  io4 = (time0 <= t4 && t4 <= time1 ? 1 : 0);

  int iolist[5] = {io0, io1, io2, io3, io4};
  double distlist[5] = {mdst0, mdst1, mdst2, mdst3, mdst4};

  // Now we only need to consider those distances for which
  // the io# flag is set to 1

  number = -1;
  mindist = -1.0;

  for (int i = 0; i < 5; i++) {
    if (iolist[i] == 1) {
      if (mindist < 0.0) {
        number = i;
        mindist = distlist[i];
      }

      if (mindist > distlist[i]) {
        number = i;
        mindist = distlist[i];
      }
    }
  }

  // Now we know the number corresponding to the boundary
  // to which the point is closest, we know the distance,
  // but we still need the normal

  distance = -1.0;
  normal_x = -1.0;
  normal_y = -1.0;

  if (number == 0) {
    distance = (dst0 < rad0 ? -mindist : mindist);

    normal_x = (x_2D - pt0x) / sqrt((x_2D - pt0x) * (x_2D - pt0x) +
                                    (y_2D - pt0y) * (y_2D - pt0y));
    normal_y = (y_2D - pt0y) / sqrt((x_2D - pt0x) * (x_2D - pt0x) +
                                    (y_2D - pt0y) * (y_2D - pt0y));
  } else if (number == 1) {
    distance = (dst1 < rad1 ? -mindist : mindist);

    normal_x = (x_2D - pt1x) / sqrt((x_2D - pt1x) * (x_2D - pt1x) +
                                    (y_2D - pt1y) * (y_2D - pt1y));
    normal_y = (y_2D - pt1y) / sqrt((x_2D - pt1x) * (x_2D - pt1x) +
                                    (y_2D - pt1y) * (y_2D - pt1y));
  } else if (number == 2) {
    distance = (dst2 < rad2 ? -mindist : mindist);

    normal_x = (x_2D - pt2x) / sqrt((x_2D - pt2x) * (x_2D - pt2x) +
                                    (y_2D - pt2y) * (y_2D - pt2y));
    normal_y = (y_2D - pt2y) / sqrt((x_2D - pt2x) * (x_2D - pt2x) +
                                    (y_2D - pt2y) * (y_2D - pt2y));
  } else if (number == 3) {
    distance = (dst3 < rad3 ? mindist : -mindist);

    normal_x = -(x_2D - pt3x) / sqrt((x_2D - pt3x) * (x_2D - pt3x) +
                                     (y_2D - pt3y) * (y_2D - pt3y));
    normal_y = -(y_2D - pt3y) / sqrt((x_2D - pt3x) * (x_2D - pt3x) +
                                     (y_2D - pt3y) * (y_2D - pt3y));
  } else if (number == 4) {
    normal_x = -1.0 / sqrt(2.0);
    normal_y = normal_x;

    if ((a - b + c - 2.0 * sqrt(2.0) * c - sqrt((b - c) * (b + 5.0 * c)) -
         x_2D - y_2D) > 0)
      distance = mindist;
    else
      distance = -mindist;
  }

  /***** Convert 2D normal to 3D coordinates *****/

  // Now that we have the normal in 2D we need to make a final
  // transformation to get it in 3D. The minimum distance stays
  // the same though. We first get the normalized direction vector.

  x = stomatocyte_3D_orientation[0];
  y = stomatocyte_3D_orientation[1];
  z = stomatocyte_3D_orientation[2];

  xd = x / sqrt(x * x + y * y + z * z);
  yd = y / sqrt(x * x + y * y + z * z);
  zd = z / sqrt(x * x + y * y + z * z);

  // We now establish the rotion matrix required to go
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

    if (zd > 0) {
      matrix[0] = 1.0;
      matrix[1] = 0.0;
      matrix[2] = 0.0;

      matrix[3] = 0.0;
      matrix[4] = 1.0;
      matrix[5] = 0.0;

      matrix[6] = 0.0;
      matrix[7] = 0.0;
      matrix[8] = 1.0;
    } else {
      matrix[0] = 1.0;
      matrix[1] = 0.0;
      matrix[2] = 0.0;

      matrix[3] = 0.0;
      matrix[4] = 1.0;
      matrix[5] = 0.0;

      matrix[6] = 0.0;
      matrix[7] = 0.0;
      matrix[8] = -1.0;
    }
  }

  // Next we determine the 3D vector between the center
  // of the stomatocyte and the point of interest

  xp = point_3D[0] - stomatocyte_3D_position[0];
  yp = point_3D[1] - stomatocyte_3D_position[1];
  zp = point_3D[2] - stomatocyte_3D_position[2];

  // Now we use the inverse matrix to find the
  // position of the point with respect to the origin
  // of the z-axis oriented stomatocyte located
  // in the origin

  xpp = matrix[0] * xp + matrix[3] * yp + matrix[6] * zp;
  ypp = matrix[1] * xp + matrix[4] * yp + matrix[7] * zp;

  // Now use this direction to orient the normal

  if (xpp * xpp + ypp * ypp > 1.e-10) {
    // The point is off the rotational symmetry
    // axis of the stomatocyte

    sin_xy = ypp / sqrt(xpp * xpp + ypp * ypp);
    cos_xy = xpp / sqrt(xpp * xpp + ypp * ypp);

    normal_x_3D = cos_xy * normal_x;
    normal_y_3D = sin_xy * normal_x;
    normal_z_3D = normal_y;
  } else {
    // The point is on the rotational symmetry
    // axis of the stomatocyte; a finite distance
    // away from the center the normal might have
    // an x and y component!

    normal_x_3D = 0.0;
    normal_y_3D = 0.0;
    normal_z_3D = (normal_y > 0.0 ? 1.0 : -1.0);
  }

  // Now we need to transform the normal back to
  // the real coordinate system

  normal_3D_x = matrix[0] * normal_x_3D + matrix[1] * normal_y_3D +
                matrix[2] * normal_z_3D;
  normal_3D_y = matrix[3] * normal_x_3D + matrix[4] * normal_y_3D +
                matrix[5] * normal_z_3D;
  normal_3D_z = matrix[6] * normal_x_3D + matrix[7] * normal_y_3D +
                matrix[8] * normal_z_3D;

  // Pass the values we obtained to ESPResSo

  if (m_direction == -1) {
    // Apply force towards inside stomatocyte

    *dist = -distance;

    vec[0] = -normal_3D_x;
    vec[1] = -normal_3D_y;
    vec[2] = -normal_3D_z;
  } else {
    // Apply force towards inside stomatocyte

    *dist = distance;

    vec[0] = normal_3D_x;
    vec[1] = normal_3D_y;
    vec[2] = normal_3D_z;
  }

  // And we are done with the stomatocyte

  vec[0] *= *dist;
  vec[1] *= *dist;
  vec[2] *= *dist;

  return 0;
}
}
