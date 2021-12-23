/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#include "config.hpp"

#ifdef LB_WALBERLA

#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "lb_interface.hpp"
#include "lb_walberla_instance.hpp"

#include <LBWalberlaNodeState.hpp>

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <functional>
#include <vector>

namespace Walberla {

void update_boundary_from_shape(std::vector<int> const &raster_flat,
                                std::vector<double> const &slip_velocity_flat) {
  lb_walberla()->update_boundary_from_shape(raster_flat, slip_velocity_flat);
  on_lb_boundary_conditions_change();
}

REGISTER_CALLBACK(update_boundary_from_shape)

void update_boundary_from_list(std::vector<int> const &nodes_flat,
                               std::vector<double> const &vel_flat) {
  lb_walberla()->update_boundary_from_list(nodes_flat, vel_flat);
  on_lb_boundary_conditions_change();
}

REGISTER_CALLBACK(update_boundary_from_list)

Utils::Vector3d get_momentum() { return lb_walberla()->get_momentum(); }

REGISTER_CALLBACK_REDUCTION(get_momentum, std::plus<>())

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos) {
  return lb_walberla()->get_velocity_at_pos(pos);
}

REGISTER_CALLBACK_ONE_RANK(get_velocity_at_pos)

boost::optional<double> get_interpolated_density_at_pos(Utils::Vector3d pos) {
  return lb_walberla()->get_interpolated_density_at_pos(pos);
}

REGISTER_CALLBACK_ONE_RANK(get_interpolated_density_at_pos)

std::size_t get_velocity_field_id() {
  return lb_walberla()->get_velocity_field_id();
}

std::size_t get_force_field_id() { return lb_walberla()->get_force_field_id(); }

} // namespace Walberla
#endif
