/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "grid_based_algorithms/lb_interpolation.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"

#include "communication.hpp"
#include "config/config.hpp"

#include <utils/Vector.hpp>

#include <iostream>
#include <stdexcept>

const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity(const Utils::Vector3d &pos) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    auto res = lb_walberla()->get_velocity_at_pos(pos / LB::get_agrid(), true);
    if (!res) {
      std::cout << this_node << ": position: [" << pos << "]\n";
      throw std::runtime_error(
          "Interpolated velocity could not be obtained from Walberla");
    }
    return *res;
#endif
  }
  throw std::runtime_error("No LB active.");
}

void lb_lbinterpolation_add_force_density(
    const Utils::Vector3d &pos, const Utils::Vector3d &force_density) {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    if (!lb_walberla()->add_force_at_pos(pos / LB::get_agrid(), force_density))
      throw std::runtime_error("Could not apply force to lb.");
#endif
  } else
    throw std::runtime_error("No LB active.");
}
