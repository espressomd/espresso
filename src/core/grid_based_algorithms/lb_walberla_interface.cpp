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

void clear_boundaries() {
  lb_walberla()->clear_boundaries();
  lb_walberla()->ghost_communication();
  on_lb_boundary_conditions_change();
}

REGISTER_CALLBACK(clear_boundaries)

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

void set_node_from_checkpoint(Utils::Vector3i ind, LBWalberlaNodeState cpt) {
  lb_walberla()->set_node_pop(ind, cpt.populations);
  lb_walberla()->set_node_last_applied_force(ind, cpt.last_applied_force);
  if (cpt.is_boundary) {
    lb_walberla()->set_node_velocity_at_boundary(ind, cpt.slip_velocity, false);
  }
}

REGISTER_CALLBACK(set_node_from_checkpoint)

boost::optional<LBWalberlaNodeState> get_node_checkpoint(Utils::Vector3i ind) {
  auto const pop = lb_walberla()->get_node_pop(ind);
  auto const laf = lb_walberla()->get_node_last_applied_force(ind);
  auto const lbb = lb_walberla()->get_node_is_boundary(ind);
  auto const vbb = lb_walberla()->get_node_velocity_at_boundary(ind);
  if (pop and laf and lbb and ((*lbb) ? vbb.has_value() : true)) {
    LBWalberlaNodeState cpnode;
    cpnode.populations = *pop;
    cpnode.last_applied_force = *laf;
    cpnode.is_boundary = *lbb;
    if (*lbb) {
      cpnode.slip_velocity = *vbb;
    }
    return cpnode;
  }
  return {boost::none};
}

REGISTER_CALLBACK_ONE_RANK(get_node_checkpoint)

void do_reallocate_ubb_field() { lb_walberla()->reallocate_ubb_field(); }

REGISTER_CALLBACK(do_reallocate_ubb_field)

void do_ghost_communication() { lb_walberla()->ghost_communication(); }

REGISTER_CALLBACK(do_ghost_communication)

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

void add_force_at_pos(Utils::Vector3d pos, Utils::Vector3d f) {
  lb_walberla()->add_force_at_pos(pos, f);
}

REGISTER_CALLBACK(add_force_at_pos)

} // namespace Walberla
#endif
