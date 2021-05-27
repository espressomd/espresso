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

#include <boost/optional.hpp>

#include <utils/Vector.hpp>

#include <vector>

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind);
boost::optional<Utils::Vector3d>
get_node_last_applied_force(Utils::Vector3i ind);
boost::optional<double> get_node_density(Utils::Vector3i ind);
boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind);
void create_vtk(unsigned delta_N, unsigned initial_count,
                unsigned flag_observables, std::string const &identifier,
                std::string const &base_folder, std::string const &prefix);
void write_vtk(std::string const &vtk_uid);
void switch_vtk(std::string const &vtk_uid, int status);
boost::optional<std::vector<double>> get_node_pop(Utils::Vector3i ind);
boost::optional<Utils::Vector6d> get_node_pressure_tensor(Utils::Vector3i ind);

void set_node_last_applied_force(Utils::Vector3i ind, Utils::Vector3d f);
void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u);
void set_ext_force_density(Utils::Vector3d f);
void set_node_density(Utils::Vector3i ind, double density);
void set_node_pop(Utils::Vector3i ind, std::vector<double> pop);
void set_node_from_checkpoint(Utils::Vector3i ind, std::vector<double> pop,
                              Utils::Vector3d f);
void do_ghost_communication();

Utils::Vector3d get_momentum();

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos);
boost::optional<double> get_interpolated_density_at_pos(Utils::Vector3d pos);

void add_force_at_pos(Utils::Vector3d pos, Utils::Vector3d f);

uint64_t get_rng_state();
void set_rng_state(uint64_t code);

} // namespace Walberla
#endif
#endif
