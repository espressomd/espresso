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
#include "lb_interface.hpp"
#include "lb_walberla_instance.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <functional>
#include <string>
#include <vector>

namespace Walberla {

boost::optional<Utils::Vector3d> get_node_velocity(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_velocity(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_velocity)

void create_vtk(unsigned delta_N, unsigned initial_count,
                unsigned flag_observables, std::string const &identifier,
                std::string const &base_folder, std::string const &prefix) {
  lb_walberla()->create_vtk(delta_N, initial_count, flag_observables,
                            identifier, base_folder, prefix);
}

REGISTER_CALLBACK(create_vtk)

void write_vtk(std::string const &vtk_uid) {
  lb_walberla()->write_vtk(vtk_uid);
}

REGISTER_CALLBACK(write_vtk)

void switch_vtk(std::string const &vtk_uid, int status) {
  lb_walberla()->switch_vtk(vtk_uid, status);
}

REGISTER_CALLBACK(switch_vtk)

boost::optional<Utils::Vector3d>
get_node_last_applied_force(Utils::Vector3i ind) {
  return lb_walberla()->get_node_last_applied_force(ind);
}

REGISTER_CALLBACK_ONE_RANK(get_node_last_applied_force)

boost::optional<double> get_node_density(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_density(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_density)

boost::optional<bool> get_node_is_boundary(Utils::Vector3i ind) {
  return lb_walberla()->get_node_is_boundary(ind);
}

REGISTER_CALLBACK_ONE_RANK(get_node_is_boundary)

boost::optional<std::vector<double>> get_node_pop(Utils::Vector3i ind) {
  return lb_walberla()->get_node_pop(ind);
}

REGISTER_CALLBACK_ONE_RANK(get_node_pop)

boost::optional<Utils::Vector6d> get_node_pressure_tensor(Utils::Vector3i ind) {
  auto res = lb_walberla()->get_node_pressure_tensor(ind);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_node_pressure_tensor)

void set_node_velocity(Utils::Vector3i ind, Utils::Vector3d u) {
  lb_walberla()->set_node_velocity(ind, u);
  lb_walberla()->ghost_communication();
}

REGISTER_CALLBACK(set_node_velocity)

void set_node_last_applied_force(Utils::Vector3i ind, Utils::Vector3d f) {
  lb_walberla()->set_node_last_applied_force(ind, f);
  lb_walberla()->ghost_communication();
}

REGISTER_CALLBACK(set_node_last_applied_force)

void set_ext_force_density(Utils::Vector3d f) {
  lb_walberla()->set_external_force(f);
}

REGISTER_CALLBACK(set_ext_force_density)

void set_node_density(Utils::Vector3i ind, double density) {
  lb_walberla()->set_node_density(ind, density);
  lb_walberla()->ghost_communication();
}

REGISTER_CALLBACK(set_node_density)

void set_node_pop(Utils::Vector3i ind, std::vector<double> pop) {
  lb_walberla()->set_node_pop(ind, pop);
  lb_walberla()->ghost_communication();
}

REGISTER_CALLBACK(set_node_pop)

void set_node_from_checkpoint(Utils::Vector3i ind, std::vector<double> pop,
                              Utils::Vector3d f) {
  lb_walberla()->set_node_pop(ind, pop);
  lb_walberla()->set_node_last_applied_force(ind, f);
}

REGISTER_CALLBACK(set_node_from_checkpoint)

void do_ghost_communication() { lb_walberla()->ghost_communication(); }

REGISTER_CALLBACK(do_ghost_communication)

Utils::Vector3d get_momentum() { return lb_walberla()->get_momentum(); }

REGISTER_CALLBACK_REDUCTION(get_momentum, std::plus<>())

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos) {
  auto res = lb_walberla()->get_velocity_at_pos(pos);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_velocity_at_pos)

boost::optional<double> get_interpolated_density_at_pos(Utils::Vector3d pos) {
  auto res = lb_walberla()->get_interpolated_density_at_pos(pos);
  return res;
}

REGISTER_CALLBACK_ONE_RANK(get_interpolated_density_at_pos)

void add_force_at_pos(Utils::Vector3d pos, Utils::Vector3d f) {
  lb_walberla()->add_force_at_pos(pos, f);
}

REGISTER_CALLBACK(add_force_at_pos)

} // namespace Walberla
#endif
