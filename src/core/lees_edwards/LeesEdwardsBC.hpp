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
#include <cmath>

struct LeesEdwardsBC {
  double pos_offset = 0.;
  double shear_velocity = 0.;
  int shear_direction = -1;
  int shear_plane_normal = -1;
  Utils::Vector3d distance(Utils::Vector3d const &d, Utils::Vector3d const &l,
                           Utils::Vector3d const &hal_l,
                           Utils::Vector3d const &l_inv,
                           std::bitset<3> const periodic) const {

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
