/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#pragma once

#include <functional>
#include <utility>

/** Pack Lees-Edwards parameters for LB. */
struct LeesEdwardsPack {
  LeesEdwardsPack(unsigned int shear_direction, unsigned int shear_plane_normal,
                  std::function<double()> get_pos_offset,
                  std::function<double()> get_shear_velocity)
      : shear_direction(shear_direction),
        shear_plane_normal(shear_plane_normal),
        get_pos_offset(std::move(get_pos_offset)),
        get_shear_velocity(std::move(get_shear_velocity)) {}
  unsigned int shear_direction;
  unsigned int shear_plane_normal;
  std::function<double()> get_pos_offset;
  std::function<double()> get_shear_velocity;
};
