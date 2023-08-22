/*
 * Copyright (C) 2019-2023 The ESPResSo project
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
#include "config/config.hpp"

#ifdef WALBERLA

#include "LBWalberla.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "system/System.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <utils/Vector.hpp>

#include <optional>

namespace LB {

double LBWalberla::get_kT() const { return lb_fluid->get_kT(); }

Utils::VectorXd<9> LBWalberla::get_pressure_tensor() const {
  return lb_fluid->get_pressure_tensor();
}

void LBWalberla::propagate() { lb_fluid->integrate(); }

void LBWalberla::lebc_sanity_checks(unsigned int shear_direction,
                                    unsigned int shear_plane_normal) const {
  lb_fluid->check_lebc(shear_direction, shear_plane_normal);
}

std::optional<Utils::Vector3d>
LBWalberla::get_velocity_at_pos(Utils::Vector3d const &pos,
                                bool consider_points_in_halo) const {
  return lb_fluid->get_velocity_at_pos(pos, consider_points_in_halo);
}

std::optional<double>
LBWalberla::get_density_at_pos(Utils::Vector3d const &pos,
                               bool consider_points_in_halo) const {
  return lb_fluid->get_density_at_pos(pos, consider_points_in_halo);
}

Utils::Vector3d LBWalberla::get_momentum() const {
  return lb_fluid->get_momentum();
}

bool LBWalberla::add_force_at_pos(Utils::Vector3d const &pos,
                                  Utils::Vector3d const &force) {
  return lb_fluid->add_force_at_pos(pos, force);
}

void LBWalberla::veto_time_step(double time_step) const {
  walberla_tau_sanity_checks("LB", lb_params->get_tau(), time_step);
}

void LBWalberla::sanity_checks() const {
  auto const agrid = lb_params->get_agrid();
  auto [my_left, my_right] = lb_fluid->get_lattice().get_local_domain();
  my_left *= agrid;
  my_right *= agrid;
  walberla_agrid_sanity_checks("LB", my_left, my_right, agrid);
  // LB time step and MD time step must agree
  walberla_tau_sanity_checks("LB", lb_params->get_tau());
}

} // namespace LB

#endif // WALBERLA
