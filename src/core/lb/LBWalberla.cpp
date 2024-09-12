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

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrate.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include "lees_edwards/lees_edwards.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <optional>
#include <variant>

namespace LB {

bool LBWalberla::is_gpu() const { return lb_fluid->is_gpu(); }

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

void LBWalberla::add_forces_at_pos(std::vector<Utils::Vector3d> const &pos,
                                   std::vector<Utils::Vector3d> const &forces) {
  lb_fluid->add_forces_at_pos(pos, forces);
}

std::vector<Utils::Vector3d>
LBWalberla::get_velocities_at_pos(std::vector<Utils::Vector3d> const &pos) {
  return lb_fluid->get_velocities_at_pos(pos);
}

void LBWalberla::veto_time_step(double time_step) const {
  walberla_tau_sanity_checks("LB", lb_params->get_tau(), time_step);
}

void LBWalberla::veto_kT(double kT) const {
  auto const energy_conversion =
      Utils::int_pow<2>(lb_params->get_agrid() / lb_params->get_tau());
  auto const lb_kT = lb_fluid->get_kT() * energy_conversion;
  if (not ::Thermostat::are_kT_equal(lb_kT, kT)) {
    throw std::runtime_error("Temperature change not supported by LB");
  }
}

void LBWalberla::sanity_checks(System::System const &system) const {
  auto const agrid = lb_params->get_agrid();
  auto [lb_left, lb_right] = lb_fluid->get_lattice().get_local_domain();
  lb_left *= agrid;
  lb_right *= agrid;
  auto const &md_left = system.local_geo->my_left();
  auto const &md_right = system.local_geo->my_right();
  walberla_agrid_sanity_checks("LB", md_left, md_right, lb_left, lb_right,
                               agrid);
  // LB time step and MD time step must agree
  walberla_tau_sanity_checks("LB", lb_params->get_tau(),
                             system.get_time_step());
}

void LBWalberla::on_lees_edwards_change() { update_collision_model(); }

void LBWalberla::update_collision_model() {
  auto const energy_conversion =
      Utils::int_pow<2>(lb_params->get_agrid() / lb_params->get_tau());
  auto const kT = lb_fluid->get_kT() * energy_conversion;
  auto const seed = lb_fluid->get_seed();
  update_collision_model(*lb_fluid, *lb_params, kT, seed);
}

void LBWalberla::update_collision_model(LBWalberlaBase &lb,
                                        LBWalberlaParams &params, double kT,
                                        unsigned int seed) {
  auto const &system = ::System::get_system();
  auto le_protocol = system.lees_edwards->get_protocol();
  if (le_protocol and
      not std::holds_alternative<LeesEdwards::Off>(*le_protocol)) {
    if (kT != 0.) {
      throw std::runtime_error(
          "Lees-Edwards LB doesn't support thermalization");
    }
    auto const &le_bc = system.box_geo->lees_edwards_bc();
    auto lees_edwards_object = std::make_unique<LeesEdwardsPack>(
        le_bc.shear_direction, le_bc.shear_plane_normal,
        [&params, le_protocol, &system]() {
          return get_pos_offset(system.get_sim_time(), *le_protocol) /
                 params.get_agrid();
        },
        [&params, le_protocol, &system]() {
          return get_shear_velocity(system.get_sim_time(), *le_protocol) *
                 (params.get_tau() / params.get_agrid());
        });
    lb.set_collision_model(std::move(lees_edwards_object));
  } else {
    lb.set_collision_model(kT, seed);
  }
  lb.ghost_communication();
}

} // namespace LB

#endif // WALBERLA
