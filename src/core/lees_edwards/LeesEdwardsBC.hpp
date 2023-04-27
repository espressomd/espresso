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
#include <initializer_list>

struct LeesEdwardsBC {
  static auto constexpr invalid_dir = 3u;
  double pos_offset = 0.;
  double shear_velocity = 0.;
  unsigned int shear_direction = invalid_dir;
  unsigned int shear_plane_normal = invalid_dir;
  Utils::Vector3d distance(Utils::Vector3d const &d, Utils::Vector3d const &l,
                           Utils::Vector3d const &hal_l,
                           Utils::Vector3d const &l_inv,
                           std::bitset<3> const periodic) const {

    Utils::Vector3d n_jumps{};
    Utils::Vector3d res = d;

    auto const le_plane_normal = static_cast<unsigned int>(shear_plane_normal);
    auto const le_direction = static_cast<unsigned int>(shear_direction);
    auto const n_le_crossings =
        std::round(res[le_plane_normal] * l_inv[le_plane_normal]);
    if (n_le_crossings >= 1.)
      res[le_direction] += pos_offset;
    if (n_le_crossings <= -1.)
      res[le_direction] -= pos_offset;

    for (auto const i : {0u, 1u, 2u}) {
      if (periodic[i]) {
        n_jumps[i] = std::round(res[i] * l_inv[i]);
        res[i] -= n_jumps[i] * l[i];
      }
    }

    return res;
  }
};

#endif
