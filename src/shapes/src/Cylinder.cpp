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

#include <cmath>
#include <shapes/Cylinder.hpp>

namespace Shapes {

std::pair<double, double> Cylinder::dist_half_pore(double r, double z) const {
  assert(z >= 0.0);
  assert(r >= 0.0);

  if (z >= m_half_length || r >= m_rad) {
    /* Outside */
    if (!m_open && z >= m_half_length && r < m_rad) {
      /* Closest feature: cap */
      return {0, -(z - m_half_length)};
    }
    if (z >= m_half_length && (m_open || r >= m_rad)) {
      /* Closest feature: ring */
      return {-(r - m_rad), -(z - m_half_length)};
    }
    /* Closest feature: cylinder */
    return {-(r - m_rad), 0};
  }
  /* Inside */
  if (!m_open && z >= m_half_length - m_rad &&
      r < (z - (m_half_length - m_rad))) {
    /* Closest feature: cap */
    return {0, m_half_length - z};
  }
  /* Closest feature: cylinder */
  return {m_rad - r, 0};
}

void Cylinder::calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const {
  /* Coordinate transform to cylinder coords
     with origin at m_center. */
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
  double dr, dz;
  std::tie(dr, dz) = dist_half_pore(r, std::abs(z));

  double side = -1;
  if (std::abs(z) >= m_half_length || r >= m_rad) /* outside */
    side = 1;

  if (z <= 0.0) {
    dz *= -1;
  }

  dist = std::sqrt(dr * dr + dz * dz) * m_direction * side;
  vec = -dr * e_r - dz * e_z;
}
} // namespace Shapes
