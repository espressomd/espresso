/*
 * Copyright (C) 2021-2022 The ESPResSo project
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
#ifndef CORE_LEES_EDWARDS_LEES_EDWARDS_BC_HPP
#define CORE_LEES_EDWARDS_LEES_EDWARDS_BC_HPP

#include <utils/Vector.hpp>

#include <bitset>
#include <cassert>
#include <cmath>

struct LeesEdwardsBC {
  double pos_offset = 0.;
  double shear_velocity = 0.;
  int shear_direction = 0;
  int shear_plane_normal = 0;
  Utils::Vector3d distance(const Utils::Vector3d &d, const Utils::Vector3d &l,
                           const Utils::Vector3d &hal_l,
                           const Utils::Vector3d &l_inv,
                           const std::bitset<3> periodic) const {
    assert(shear_plane_normal != shear_direction);
    assert(shear_direction >= 0 and shear_direction <= 2);
    assert(shear_plane_normal >= 0 and shear_plane_normal <= 2);

    Utils::Vector3d n_jumps{};
    Utils::Vector3d res = d;

    double n_le_crossings =
        std::round(res[shear_plane_normal] * l_inv[shear_plane_normal]);
    res[shear_direction] += n_le_crossings * pos_offset;

    for (int i : {0, 1, 2}) {
      if (periodic[i])
        n_jumps[i] = std::round(res[i] * l_inv[i]);
      res[i] -= n_jumps[i] * l[i];
    }

    return res;
  }
};

#endif
