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

#include <shapes/SimplePore.hpp>

#include <cassert>

namespace Shapes {
/**
 * @brief Calculate the distance function in the coordinates of the cylinder.
 *
 * @param r Distance from the cylinder axis.
 * @param z Distance from the center along the axis.
 *
 * @returns The distance vector from the surface in the cylinder system.
 */
std::pair<double, double> SimplePore::dist_half_pore(double r, double z) const {
  assert(z >= 0.0);
  assert(r >= 0.0);

  /*
   *  We have to find the line that splits area 1 (r determines distance) from
   *  area 2 (z determines distance) inside pore. In area 3 we have to consider
   *  z and r to determine the distance.
   *
   *   |        x
   *   |   2  x
   *   |    x
   *  _|_ x    1
   *    \|       ^ r
   *  3  |-------|
   *     |   z <-
   */

  if ((z <= c_z) && (r <= (c_z + c_r - z))) {
    /* Cylinder section, inner */
    return {m_rad - r, 0};
  }
  if (((z >= c_z) && (r >= c_r)) || ((z <= c_z) && (r > (c_z + c_r - z)))) {
    /* Wall section and outer cylinder */
    return {0, m_half_length - z};
  }
  /* Smoothing area */
  /* Vector to center of torus segment */
  auto const dr = c_r - r;
  auto const dz = c_z - z;

  /* Rescale to surface */
  auto const d = std::sqrt(dr * dr + dz * dz);
  auto const fac = (d - m_smoothing_rad) / d;

  return {fac * dr, fac * dz};
}

void SimplePore::calculate_dist(const Utils::Vector3d &pos, double &dist,
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
  if (((dz == 0) && (r <= m_rad)) ||                  // cylinder section
      ((dr == 0) && (std::abs(z) > m_half_length))) { // ||
    side = 1;
  } else {
    // smoothing area
    if (std::abs(z) >= c_z) {
      auto const angle = std::asin((std::abs(z) - c_z) / m_smoothing_rad);
      auto const dist_offset =
          m_smoothing_rad - (std::cos(angle) * m_smoothing_rad);
      if (m_half_length < std::abs(z) || r <= (m_rad + dist_offset)) {
        side = 1;
      }
    }
  }

  if (z <= 0.0) {
    dz *= -1;
  }

  dist = std::sqrt(dr * dr + dz * dz) * side;
  vec = -dr * e_r - dz * e_z;
}
} // namespace Shapes
