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

#include "HollowCone.hpp"

#include <cmath>
#include <cassert>

using namespace std;

namespace Shapes {
int HollowCone::calculate_dist(const double *ppos, double *dist, double *vec) const {
  int number = -1;
  double r0, r1, w, alpha, xd, yd, zd, mu, x_2D, y_2D, t0, t1, t2, time1, time2,
      time3, time4, mdst0, mdst1, mdst2, mdst3, mindist, normalize, x, y, z,
      distance, normal_x, normal_y, direction = 0.0, xp, yp, zp, xpp, ypp, sin_xy,
      cos_xy, normal_x_3D, normal_y_3D, normal_z_3D, normal_3D_x, normal_3D_y,
      normal_3D_z;

  double closest_point_3D[3] = {-1.0, -1.0, -1.0};

  // Set the dimensions of the hollow cone

  r0 = m_inner_radius;
  r1 = m_outer_radius;
  w = m_width;
  alpha = m_opening_angle;

  // Set the position and orientation of the hollow cone

  double hollow_cone_3D_position[3] = {m_position_x, m_position_y, m_position_z};

  double hollow_cone_3D_orientation[3] = {m_orientation_x, m_orientation_y,
                                          m_orientation_z};

  // Set the point for which we want to know the distance

  double point_3D[3];

  point_3D[0] = ppos[0];
  point_3D[1] = ppos[1];
  point_3D[2] = ppos[2];

  /***** Convert 3D coordinates to 2D planar coordinates *****/

  // Calculate the point on position + mu * orientation,
  // where the difference segment is orthogonal

  mu = (hollow_cone_3D_orientation[0] * point_3D[0] +
        hollow_cone_3D_orientation[1] * point_3D[1] +
        hollow_cone_3D_orientation[2] * point_3D[2] -
        hollow_cone_3D_position[0] * hollow_cone_3D_orientation[0] -
        hollow_cone_3D_position[1] * hollow_cone_3D_orientation[1] -
        hollow_cone_3D_position[2] * hollow_cone_3D_orientation[2]) /
       (hollow_cone_3D_orientation[0] * hollow_cone_3D_orientation[0] +
        hollow_cone_3D_orientation[1] * hollow_cone_3D_orientation[1] +
        hollow_cone_3D_orientation[2] * hollow_cone_3D_orientation[2]);

  // Then the closest point to the line is

  closest_point_3D[0] =
      hollow_cone_3D_position[0] + mu * hollow_cone_3D_orientation[0];
  closest_point_3D[1] =
      hollow_cone_3D_position[1] + mu * hollow_cone_3D_orientation[1];
  closest_point_3D[2] =
      hollow_cone_3D_position[2] + mu * hollow_cone_3D_orientation[2];

  // So the shortest distance to the line is

  x_2D = sqrt((closest_point_3D[0] - point_3D[0]) *
                  (closest_point_3D[0] - point_3D[0]) +
              (closest_point_3D[1] - point_3D[1]) *
                  (closest_point_3D[1] - point_3D[1]) +
              (closest_point_3D[2] - point_3D[2]) *
                  (closest_point_3D[2] - point_3D[2]));

  y_2D =
      mu * sqrt(hollow_cone_3D_orientation[0] * hollow_cone_3D_orientation[0] +
                hollow_cone_3D_orientation[1] * hollow_cone_3D_orientation[1] +
                hollow_cone_3D_orientation[2] * hollow_cone_3D_orientation[2]);

  /***** Use the obtained planar coordinates in distance function *****/

  // Calculate intermediate results which we need to determine
  // the distance

  t0 = (y_2D * cos(alpha) + (x_2D - r0) * sin(alpha)) / r1;
  t1 = (w - 2.0 * (x_2D - r0) * cos(alpha) + 2.0 * y_2D * sin(alpha)) /
       (2.0 * w);
  t2 = (w + 2.0 * (x_2D - r0) * cos(alpha) - 2.0 * y_2D * sin(alpha)) /
       (2.0 * w);

  if (t0 >= 0.0 && t0 <= 1.0) {
    time1 = t0;
    time2 = t0;
  } else if (t0 > 1.0) {
    time1 = 1.0;
    time2 = 1.0;
  } else {
    time1 = 0.0;
    time2 = 0.0;
  }

  if (t1 >= 0.0 && t1 <= 1.0) {
    time3 = t1;
  } else if (t1 > 1.0) {
    time3 = 1.0;
  } else {
    time3 = 0.0;
  }

  if (t2 >= 0.0 && t2 <= 1.0) {
    time4 = t2;
  } else if (t2 > 1.0) {
    time4 = 1.0;
  } else {
    time4 = 0.0;
  }

  mdst0 =
      x_2D * x_2D + y_2D * y_2D - 2.0 * x_2D * r0 + r0 * r0 +
      r1 * r1 * time1 * time1 + 0.25 * w * w +
      (-2.0 * y_2D * r1 * time1 + (x_2D - r0) * w) * cos(alpha) -
      (2.0 * x_2D * r1 * time1 - 2.0 * r0 * r1 * time1 + y_2D * w) * sin(alpha);
  mdst0 = sqrt(mdst0);

  mdst1 = x_2D * x_2D + y_2D * y_2D - 2.0 * x_2D * r0 + r0 * r0 +
          r1 * r1 * time2 * time2 + 0.25 * w * w +
          (-2.0 * y_2D * r1 * time2 + (-x_2D + r0) * w) * cos(alpha) +
          (-2.0 * x_2D * r1 * time2 + 2.0 * r0 * r1 * time2 + y_2D * w) *
              sin(alpha);
  mdst1 = sqrt(mdst1);

  mdst2 = x_2D * x_2D + y_2D * y_2D - 2.0 * x_2D * r0 + r0 * r0 + 0.25 * w * w -
          time3 * w * w + time3 * time3 * w * w +
          (x_2D - r0) * (-1.0 + 2.0 * time3) * w * cos(alpha) -
          y_2D * (-1.0 + 2.0 * time3) * w * sin(alpha);
  mdst2 = sqrt(mdst2);

  mdst3 = x_2D * x_2D + y_2D * y_2D - 2.0 * x_2D * r0 + r0 * r0 + 0.25 * w * w -
          time4 * w * w + time4 * time4 * w * w -
          (x_2D - r0) * (-1.0 + 2.0 * time4) * w * cos(alpha) +
          y_2D * (-1.0 + 2.0 * time4) * w * sin(alpha) + r1 * r1 -
          2.0 * y_2D * r1 * cos(alpha) +
          (-2.0 * x_2D * r1 + 2.0 * r0 * r1) * sin(alpha);
  mdst3 = sqrt(mdst3);

  double distlist[4] = {mdst0, mdst1, mdst2, mdst3};

  // Now we only need to determine which distance is minimal
  // and remember which one it is

  mindist = -1.0;

  for (int i = 0; i < 4; i++) {
    if (mindist < 0.0) {
      number = i;
      mindist = distlist[i];
    }

    if (mindist > distlist[i]) {
      number = i;
      mindist = distlist[i];
    }
  }

  // Now we know the number corresponding to the boundary
  // to which the point is closest, we know the distance,
  // but we still need the normal

  distance = -1.0;
  normal_x = -1.0;
  normal_y = -1.0;

  if (number == 0) {
    normal_x = x_2D - r0 + 0.5 * w * cos(alpha) - r1 * time1 * sin(alpha);
    normal_y = y_2D - r1 * time1 * cos(alpha) - 0.5 * w * sin(alpha);
    normalize = 1.0 / sqrt(normal_x * normal_x + normal_y * normal_y);
    normal_x *= normalize;
    normal_y *= normalize;

    direction = -normal_x * cos(alpha) + normal_y * sin(alpha);

    if (fabs(direction) < 1.0e-06 &&
        (fabs(time1 - 0.0) < 1.0e-06 || fabs(time1 - 1.0) < 1.0e-06)) {
      if (fabs(time1 - 0.0) < 1.0e-06) {
        direction = -normal_x * sin(alpha) - normal_y * cos(alpha);
      } else {
        direction = normal_x * sin(alpha) + normal_y * cos(alpha);
      }
    }

    if (direction > 0.0) {
      distance = mindist;
    } else {
      distance = -mindist;
      normal_x *= -1.0;
      normal_y *= -1.0;
    }
  } else if (number == 1) {
    normal_x = x_2D - r0 - 0.5 * w * cos(alpha) - r1 * time2 * sin(alpha);
    normal_y = y_2D - r1 * time2 * cos(alpha) + 0.5 * w * sin(alpha);
    normalize = 1.0 / sqrt(normal_x * normal_x + normal_y * normal_y);
    normal_x *= normalize;
    normal_y *= normalize;

    direction = normal_x * cos(alpha) - normal_y * sin(alpha);

    if (fabs(direction) < 1.0e-06 &&
        (fabs(time2 - 0.0) < 1.0e-06 || fabs(time2 - 1.0) < 1.0e-06)) {
      if (fabs(time2 - 0.0) < 1.0e-06) {
        direction = -normal_x * sin(alpha) - normal_y * cos(alpha);
      } else {
        direction = normal_x * sin(alpha) + normal_y * cos(alpha);
      }
    }

    if (direction > 0.0) {
      distance = mindist;
    } else {
      distance = -mindist;
      normal_x *= -1.0;
      normal_y *= -1.0;
    }
  } else if (number == 2) {
    normal_x = x_2D - r0 - 0.5 * (1.0 - 2.0 * time3) * w * cos(alpha);
    normal_y = y_2D + 0.5 * (1.0 - 2.0 * time3) * w * sin(alpha);
    normalize = 1.0 / sqrt(normal_x * normal_x + normal_y * normal_y);
    normal_x *= normalize;
    normal_y *= normalize;

    direction = -normal_x * sin(alpha) - normal_y * cos(alpha);

    if (fabs(direction) < 1.0e-06 &&
        (fabs(time3 - 0.0) < 1.0e-06 || fabs(time3 - 1.0) < 1.0e-06)) {
      if (fabs(time3 - 0.0) < 1.0e-06) {
        direction = normal_x * cos(alpha) - normal_y * sin(alpha);
      } else {
        direction = -normal_x * cos(alpha) + normal_y * sin(alpha);
      }
    }

    if (direction > 0.0) {
      distance = mindist;
    } else {
      distance = -mindist;
      normal_x *= -1.0;
      normal_y *= -1.0;
    }
  } else if (number == 3) {
    normal_x = x_2D - r0 + 0.5 * (1.0 - 2.0 * time4) * w * cos(alpha) -
               r1 * sin(alpha);
    normal_y =
        y_2D - 0.5 * (1.0 - 2.0 * time4) * w * sin(alpha) - r1 * cos(alpha);
    normalize = 1.0 / sqrt(normal_x * normal_x + normal_y * normal_y);
    normal_x *= normalize;
    normal_y *= normalize;

    direction = normal_x * sin(alpha) + normal_y * cos(alpha);

    if (fabs(direction) < 1.0e-06 &&
        (fabs(time4 - 0.0) < 1.0e-06 || fabs(time4 - 1.0) < 1.0e-06)) {
      if (fabs(time4 - 0.0) < 1.0e-06) {
        direction = -normal_x * cos(alpha) + normal_y * sin(alpha);
      } else {
        direction = normal_x * cos(alpha) - normal_y * sin(alpha);
      }
    }

    if (direction > 0.0) {
      distance = mindist;
    } else {
      distance = -mindist;
      normal_x *= -1.0;
      normal_y *= -1.0;
    }
  } else {
    assert(0);
  }

  /***** Convert 2D normal to 3D coordinates *****/

  // Now that we have the normal in 2D we need to make a final
  // transformation to get it in 3D. The minimum distance stays
  // the same though. We first get the normalized direction vector.

  x = hollow_cone_3D_orientation[0];
  y = hollow_cone_3D_orientation[1];
  z = hollow_cone_3D_orientation[2];

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
  // of the hollow cylinder and the point of interest

  xp = point_3D[0] - hollow_cone_3D_position[0];
  yp = point_3D[1] - hollow_cone_3D_position[1];
  zp = point_3D[2] - hollow_cone_3D_position[2];

  // Now we use the inverse matrix to find the
  // position of the point with respect to the origin
  // of the z-axis oriented hollow cone located
  // in the origin

  xpp = matrix[0] * xp + matrix[3] * yp + matrix[6] * zp;
  ypp = matrix[1] * xp + matrix[4] * yp + matrix[7] * zp;

  // Now use this direction to orient the normal

  if (xpp * xpp + ypp * ypp > 1.0e-10) {
    // The point is off the rotational symmetry
    // axis of the hollow cone

    sin_xy = ypp / sqrt(xpp * xpp + ypp * ypp);
    cos_xy = xpp / sqrt(xpp * xpp + ypp * ypp);

    normal_x_3D = cos_xy * normal_x;
    normal_y_3D = sin_xy * normal_x;
    normal_z_3D = normal_y;
  } else {
    // The point is on the rotational symmetry
    // axis of the hollow cone; a finite distance
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

  if (direction == -1) {
    // Apply force towards inside hollow cone

    *dist = -distance;

    vec[0] = -normal_3D_x;
    vec[1] = -normal_3D_y;
    vec[2] = -normal_3D_z;
  } else {
    // Apply force towards inside hollow cone

    *dist = distance;

    vec[0] = normal_3D_x;
    vec[1] = normal_3D_y;
    vec[2] = normal_3D_z;
  }

  // And we are done with the hollow cone

  vec[0] *= *dist;
  vec[1] *= *dist;
  vec[2] *= *dist;

  return 0;
}
}
