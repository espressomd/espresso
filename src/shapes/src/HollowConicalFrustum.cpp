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
#include <boost/algorithm/clamp.hpp>

#include <shapes/HollowConicalFrustum.hpp>
#include <utils/Vector.hpp>
#include <utils/math/coordinate_transformation.hpp>

namespace Shapes {

void HollowConicalFrustum::calculate_dist(const Utils::Vector3d &pos,
                                          double &dist,
                                          Utils::Vector3d &vec) const {
  // transform given position to cylindrical coordinates in the reference frame
  // of the cone
  auto const pos_cyl =
      Utils::transform_coordinate_cartesian_to_cylinder(pos - m_center, m_axis);
  // clang-format off
  /*
   * the following implementation is based on:
   *   - defining the cone in the cylindrical 2d coordinates with e_z = m_axis
   * and r in the plane of axis and pos (1)
   *   - defining the normal to the cone (2)
   *   - find the intersection between (1) and (2)
   *
   *   r_cone(z) = m * z + (r1 + r2) / 2, -l/2 <= z <= l/2 with m = (r1 - r2) / l
   *   r_normal(z) = -1/m * z + r_pos + 1 / m * z_pos
   *   r_cone = r_normal => z_intersection = (-m*(r1+r2-2*r_pos)+2*z_pos) / (2*(1+m**2))
   */
  // clang-format on
  auto const m = (m_r1 - m_r2) / m_length;
  auto z_intersection = (-m * (m_r1 + m_r2 - 2 * pos_cyl[0]) + 2 * pos_cyl[2]) /
                        (2 * (1 + m * m));
  // Limit the possible z values to the range in which the cone exists.
  z_intersection =
      boost::algorithm::clamp(z_intersection, -m_length / 2.0, m_length / 2.0);
  auto r_cone = [this](auto z) {
    return (m_r1 - m_r2) / m_length * z + (m_r1 + m_r2) / 2.0;
  };
  auto const r_intersection = r_cone(z_intersection);
  // Transform back to cartesian coordinates.
  auto const pos_intersection =
      Utils::transform_coordinate_cylinder_to_cartesian(
          {r_intersection, pos_cyl[1], z_intersection}, m_axis) +
      m_center;

  auto const u = (pos - pos_intersection).normalize();
  auto const d = (pos - pos_intersection).norm() - 0.5 * m_thickness;
  auto const in_out = (d < 0.0) ? -m_direction : m_direction;
  dist = in_out * std::abs(d);
  vec = d * u;
}
} // namespace Shapes
