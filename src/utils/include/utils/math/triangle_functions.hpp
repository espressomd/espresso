/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef UTILS_MATH_TRIANGEL_FUNCTIONS_HPP
#define UTILS_MATH_TRIANGEL_FUNCTIONS_HPP

#include "utils/Vector.hpp"

#include <boost/algorithm/clamp.hpp>

#include <cmath>

namespace Utils {
/**
 * @brief Computes the normal vector of a triangle.
 *
 * The sign convention is such that P1P2, P1P3 and
 * the normal form a right-handed system.
 * The normal vector is not normalized, e.g. its length
 * is arbitrary.
 */
inline Vector3d get_n_triangle(const Vector3d &P1, const Vector3d &P2,
                               const Vector3d &P3) {
  auto const u = P2 - P1;
  auto const v = P3 - P1;

  return vector_product(u, v);
}

/** Computes the area of triangle between vectors P1,P2,P3, by computing
 *  the cross product P1P2 x P1P3 and taking the half of its norm.
 */
inline double area_triangle(const Vector3d &P1, const Vector3d &P2,
                            const Vector3d &P3) {
  return 0.5 * get_n_triangle(P1, P2, P3).norm();
}

/** @brief Returns the angle between two triangles in 3D space


Returns the angle between two triangles in 3D space given by points P1P2P3 and
P2P3P4. Note, that the common edge is given as the second and the third
argument. Here, the angle can have values from 0 to 2 * PI, depending on the
orientation of the two triangles. So the angle can be convex or concave. For
each triangle, an inward direction has been defined as the direction of one of
the two normal vectors. Particularly, for triangle P1P2P3 it is the vector N1 =
P2P1 x P2P3 and for triangle P2P3P4 it is N2 = P2P3 x P2P4. The method first
computes the angle between N1 and N2, which gives always value between 0 and PI
and then it checks whether this value must be corrected to a value between PI
and 2 * PI.

As an example, consider 4 points A,B,C,D in space given by coordinates A =
[1,1,1], B = [2,1,1], C = [1,2,1], D = [1,1,2]. We want to determine the angle
between triangles ABC and ACD. In case that the orientations of the triangle ABC
is [0,0,1] and orientation of ACD is [1,0,0], then the resulting angle must be
PI/2.0. To get correct result, note that the common edge is AC, and one must
call the method as angle_btw_triangles(B,A,C,D). With this call we have ensured
that N1 = AB x AC (which coincides with [0,0,1]) and N2 = AC x AD (which
coincides with [1,0,0]). Alternatively, if the orientations of the two triangles
were the opposite, the correct call would be angle_btw_triangles(B,C,A,D) so
that N1 = CB x CA and N2 = CA x CD.

*/
inline double angle_btw_triangles(const Vector3d &P1, const Vector3d &P2,
                                  const Vector3d &P3, const Vector3d &P4) {
  auto const normal1 = get_n_triangle(P2, P1, P3);
  auto const normal2 = get_n_triangle(P2, P3, P4);
  auto const cosine = boost::algorithm::clamp(
      normal1 * normal2 / std::sqrt(normal1.norm2() * normal2.norm2()), -1.0,
      1.0);
  // The angle between the faces (not considering
  // the orientation, always less or equal to Pi)
  // is equal to Pi minus angle between the normals
  auto const phi = M_PI - std::acos(cosine);

  // Now we need to determine, if the angle between two triangles is less than
  // Pi or more than Pi. To do this, we check
  // if the point P4 lies in the halfspace given by triangle P1P2P3 and the
  // normal to this triangle. If yes, we have
  // angle less than Pi, if not, we have angle more than Pi.
  // General equation of the plane is n_x*x + n_y*y + n_z*z + d = 0 where
  // (n_x,n_y,n_z) is the normal to the plane.
  // Point P1 lies in the plane, therefore d = -(n_x*P1_x + n_y*P1_y + n_z*P1_z)
  // Point P4 lies in the halfspace given by normal iff n_x*P4_x + n_y*P4_y +
  // n_z*P4_z + d >= 0
  if (normal1 * P4 - normal1 * P1 < 0)
    return 2 * M_PI - phi;

  return phi;
}
} // namespace Utils

#endif
