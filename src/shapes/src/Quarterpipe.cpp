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

#include <shapes/Quarterpipe.hpp>
#include <utils/Vector.hpp>
#include <utils/math/abs.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>

namespace Shapes {

void Quarterpipe::calculate_dist(const Utils::Vector3d &pos, double &dist,
                                 Utils::Vector3d &vec) const {

  // Please have a look at the Doxygen page of this shape for a detailed sketch
  // and explanation of terms used below

  auto pos_shifted = pos - m_center;
  auto rel_z = pos_shifted * m_axis;

  // first, do only the 2d problem, add the axis projection later

  auto pos_projected = pos_shifted - m_axis * rel_z;
  auto pos_norm = pos_projected.norm();

  if (pos_norm == 0.) {
    dist = m_radius;
    vec = -dist * m_orientation;
  } else {
    // go to cylindrical coordinates (r,delta_phi) where m_axis is the z-axis
    // and  m_orientation defines the delta_phi = 0 axis
    auto x = pos_projected * m_orientation;
    auto y = pos_projected * Utils::vector_product(m_axis, m_orientation);
    auto pos_delta_phi = std::atan2(y, x);

    // point 1 and 2 are the end points of the shape where the cylinder
    // intersects the cuboid. Their position vectors also define normals of the
    // cuboid side planes 1 and 2
    auto to_point_1 =
        m_radius * vec_rotate(m_axis, -0.5 * m_delta_phi, m_orientation);
    auto to_point_2 =
        m_radius * vec_rotate(m_axis, 0.5 * m_delta_phi, m_orientation);

    if (pos_delta_phi <= -0.5 * m_delta_phi) {
      // Zone 1
      vec = pos_projected - to_point_1;
      dist = vec.norm();
    } else if (pos_delta_phi >= 0.5 * m_delta_phi) {
      // Zone 2
      vec = pos_projected - to_point_2;
      dist = vec.norm();
    } else if (pos_norm < m_radius) {
      // Zone 3
      dist = m_radius - pos_norm;
      vec = -dist * pos_projected / pos_norm;
    } else {
      // Zones 4 to 9 are hard to distingish, so all distances are calculated
      // and the minimum chosen

      // Zone 4, 5, 7 or 9
      auto dist_plane_1 =
          (pos_projected * to_point_1 / to_point_1.norm()) - to_point_1.norm();
      auto vec_plane_1 = dist_plane_1 * to_point_1 / to_point_1.norm();

      // Zone 4, 6, 8 or 9
      auto dist_plane_2 =
          (pos_projected * to_point_2 / to_point_2.norm()) - to_point_2.norm();
      auto vec_plane_2 = dist_plane_2 * to_point_2 / to_point_2.norm();

      // Distinguish the Zones
      if (dist_plane_1 >= 0. and dist_plane_2 >= 0.) {
        // Zone 9
        auto corner_point = std::sqrt(2.) * m_radius * m_orientation;
        vec = pos_projected - corner_point;
        dist = vec.norm();
      } else if (dist_plane_1 < 0. and dist_plane_2 < 0.) {
        // Zone 4, 5 or 6

        // First, find which side is closer
        // Distinguishes zones 4 and 5 from 4 and 6
        bool in_zone_4_or_5 = dist_plane_1 > dist_plane_2;
        auto dist_plane = in_zone_4_or_5 ? dist_plane_1 : dist_plane_2;
        auto vec_plane = in_zone_4_or_5 ? vec_plane_1 : vec_plane_2;

        // Now find if in 4 or 5/6
        auto dist_cyl = m_radius - pos_norm;
        auto vec_cyl = -dist_cyl * pos_projected / pos_norm;

        // both are negative
        bool in_zone_4 = dist_cyl > dist_plane;
        dist = in_zone_4 ? dist_cyl : dist_plane;
        vec = in_zone_4 ? vec_cyl : vec_plane;
      } else {
        // Zone 7 or 8
        bool in_zone_7 = dist_plane_1 >= 0;
        dist = in_zone_7 ? dist_plane_1 : dist_plane_2;
        vec = in_zone_7 ? vec_plane_1 : vec_plane_2;
      }
    }
  } // 2D end, now to the 3rd dimension

  auto upper_lower = rel_z > 0 ? 1 : -1;
  // signed distance to the closer top/bottom wall
  auto dist_cover = upper_lower * rel_z - 0.5 * m_height;

  // If the point is in zone A or B, nothing has to change

  if (dist_cover < 0. and dist < dist_cover) {
    // Zone C (abs(dist_cover) < abs(dist))
    dist = dist_cover;
    vec = upper_lower * dist * m_axis;
  } else if (dist_cover > 0) {
    // Zone D or E
    auto z_to_shape = rel_z - upper_lower * 0.5 * m_height;
    if (dist < 0) {
      // Zone D
      dist = Utils::abs(z_to_shape);
      vec = upper_lower * dist * m_axis;
    } else {
      // Zone E
      vec += z_to_shape * m_axis;
      dist = vec.norm();
    }
  }
}
} // namespace Shapes
