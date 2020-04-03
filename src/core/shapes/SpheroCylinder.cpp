/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "SpheroCylinder.hpp"
#include <cmath>

namespace Shapes {

void SpheroCylinder::calculate_dist(const Utils::Vector3d &pos, double &dist,
                                    Utils::Vector3d &vec) const {
  /* Coordinate transform to cylinder coords
     with origin at m_center. */

  /* Vector cylinder center<->particle */
  Utils::Vector3d const c_dist = pos - m_center;

  auto const z = e_z * c_dist;
  auto const r_vec = c_dist - z * e_z;
  auto const r = r_vec.norm();

  /* If exactly on the axis, chose e_r orthogonal
     to e_z. */
  auto const e_r = (r == 0) ? e_r_axis : r_vec / r;

  /* The pore has mirror symmetry in z with regard to
     the center in the {r,z} system. We calculate always
     for the z > 0 case, and flip the result if appropriate. */
  double dr;
  double z_abs = std::abs(z);
  double side = 1;

  if (r >= m_rad ||
      (z_abs >= m_half_length &&
       std::sqrt(r * r + std::pow(z_abs - m_half_length, 2)) > m_rad)) {
    /* Outside */
    if (z_abs >= m_half_length) {
      /* Closest feature: hemisphere */
      double dir = 1;
      if (z < 0)
        dir = -1;
      Utils::Vector3d c_dist_cap = pos - (m_center + dir * e_z * m_half_length);
      dist = c_dist_cap.norm() - m_rad;
      c_dist_cap.normalize();
      vec = dist * c_dist_cap;
      dist *= m_direction;
      return;
    }
    /* Closest feature: cylinder */
    dr = -(r - m_rad);

  } else {
    side = -1;
    /* Inside */
    if (z_abs <= m_half_length) {
      /* Closest feature: cylinder */
      dr = m_rad - r;
    } else {
      /* Closest feature: hemisphere */
      double dir = 1;
      if (z < 0)
        dir = -1;
      Utils::Vector3d c_dist_cap =
          -(pos - (m_center + dir * e_z * m_half_length));
      dist = m_rad - c_dist_cap.norm();
      c_dist_cap.normalize();
      vec = dist * c_dist_cap;
      dist *= -m_direction;
      return;
    }
  }

  dist = std::abs(dr) * m_direction * side;
  vec = -dr * e_r;
}
} // namespace Shapes
