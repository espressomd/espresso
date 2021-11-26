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
#ifndef LB_WALBERLA_INTERFACE_HPP
#define LB_WALBERLA_INTERFACE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <vector>

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind);
boost::optional<Utils::Vector3d>
get_node_velocity_at_boundary(Utils::Vector3i ind);
boost::optional<Utils::Vector3d>
get_node_last_applied_force(Utils::Vector3i ind);
boost::optional<double> get_node_density(Utils::Vector3i ind);
boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind);
void clear_boundaries();
void update_boundary_from_shape(std::vector<int> const &raster_flat,
                                std::vector<double> const &slip_velocity_flat);
void update_boundary_from_list(std::vector<int> const &nodes_flat,
                               std::vector<double> const &vel_flat);
boost::optional<Utils::Vector3d> get_node_boundary_force(Utils::Vector3i ind);
void remove_node_from_boundary(Utils::Vector3i ind);
boost::optional<std::vector<double>> get_node_pop(Utils::Vector3i ind);
boost::optional<Utils::Vector6d> get_node_pressure_tensor(Utils::Vector3i ind);

void set_node_last_applied_force(Utils::Vector3i ind, Utils::Vector3d f);
void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u);
void set_node_velocity_at_boundary(Utils::Vector3i ind, Utils::Vector3d u);
void set_node_density(Utils::Vector3i ind, double density);
void set_node_pop(Utils::Vector3i ind, std::vector<double> pop);
void set_node_from_checkpoint(Utils::Vector3i ind, std::vector<double> pop,
                              Utils::Vector3d f, Utils::Vector3d v, bool lbb);
void do_reallocate_ubb_field();
void do_ghost_communication();

Utils::Vector3d get_momentum();

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos);
boost::optional<double> get_interpolated_density_at_pos(Utils::Vector3d pos);

void add_force_at_pos(Utils::Vector3d pos, Utils::Vector3d f);

} // namespace Walberla
#endif
#endif
