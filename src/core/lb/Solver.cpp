/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "lb/Implementation.hpp"
#include "lb/Solver.hpp"
#include "lb/utils.hpp"

#include "lb/LBNone.hpp"
#include "lb/LBWalberla.hpp"

#include "BoxGeometry.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#ifdef WALBERLA
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#endif

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace LB {

Solver::Solver() { impl = std::make_unique<Implementation>(); }

static auto is_solver_set(std::unique_ptr<Solver::Implementation> const &ptr) {
  return ptr != nullptr and ptr->solver.has_value();
}

static void check_solver(std::unique_ptr<Solver::Implementation> const &ptr) {
  if (not is_solver_set(ptr)) {
    throw NoLBActive{};
  }
}

bool Solver::is_solver_set() const { return LB::is_solver_set(impl); }

void Solver::reset() {
  System::get_system().lb.impl->solver = std::nullopt;
  m_conv = Conversions{};
}

void Solver::propagate() {
  check_solver(impl);
  std::visit([](auto &ptr) { ptr->propagate(); }, *impl->solver);
}

void Solver::sanity_checks() const {
  if (impl->solver) {
    auto const &system = get_system();
    std::visit([&](auto &ptr) { ptr->sanity_checks(system); }, *impl->solver);
  }
}

void Solver::veto_time_step(double time_step) const {
  if (impl->solver) {
    std::visit([=](auto &ptr) { ptr->veto_time_step(time_step); },
               *impl->solver);
  }
}

void Solver::veto_kT(double kT) const {
  if (impl->solver) {
    std::visit([=](auto &ptr) { ptr->veto_kT(kT); }, *impl->solver);
  }
}

void Solver::lebc_sanity_checks(unsigned int shear_direction,
                                unsigned int shear_plane_normal) const {
  if (impl->solver) {
    auto const callback = [=](auto &ptr) {
      ptr->lebc_sanity_checks(shear_direction, shear_plane_normal);
    };
    std::visit(callback, *impl->solver);
  }
}

void Solver::on_cell_structure_change() {
  if (impl->solver) {
    auto &solver = *impl->solver;
    std::visit([](auto &ptr) { ptr->on_cell_structure_change(); }, solver);
  }
}

void Solver::veto_boxl_change() const {
  if (impl->solver) {
    std::visit([](auto const &ptr) { ptr->veto_boxl_change(); }, *impl->solver);
  }
}

void Solver::on_boxl_change() {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->on_boxl_change(); }, *impl->solver);
  }
}

void Solver::on_node_grid_change() {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->on_node_grid_change(); }, *impl->solver);
  }
}

void Solver::on_timestep_change() {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->on_timestep_change(); }, *impl->solver);
  }
}

void Solver::on_temperature_change() {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->on_temperature_change(); }, *impl->solver);
  }
}

bool Solver::is_gpu() const {
  check_solver(impl);
  return std::visit([](auto &ptr) { return ptr->is_gpu(); }, *impl->solver);
}

double Solver::get_agrid() const {
  check_solver(impl);
  return std::visit([](auto &ptr) { return ptr->get_agrid(); }, *impl->solver);
}

double Solver::get_tau() const {
  check_solver(impl);
  return std::visit([](auto &ptr) { return ptr->get_tau(); }, *impl->solver);
}

double Solver::get_kT() const {
  check_solver(impl);
  return std::visit([](auto &ptr) { return ptr->get_kT(); }, *impl->solver);
}

Utils::VectorXd<9> Solver::get_pressure_tensor() const {
  check_solver(impl);
  return std::visit([](auto &ptr) { return ptr->get_pressure_tensor(); },
                    *impl->solver);
}

std::optional<Utils::Vector3d>
Solver::get_interpolated_velocity(Utils::Vector3d const &pos) const {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  return std::visit(
      [&](auto &ptr) {
        auto const &box_geo = *System::get_system().box_geo;
        auto const lb_pos = box_geo.folded_position(pos) * m_conv.pos_to_lb;
        return ptr->get_velocity_at_pos(lb_pos, false);
      },
      *impl->solver);
}

std::optional<double>
Solver::get_interpolated_density(Utils::Vector3d const &pos) const {
  return std::visit(
      [&](auto &ptr) {
        auto const &box_geo = *System::get_system().box_geo;
        auto const lb_pos = box_geo.folded_position(pos) * m_conv.pos_to_lb;
        return ptr->get_density_at_pos(lb_pos, false);
      },
      *impl->solver);
}

Utils::Vector3d
Solver::get_coupling_interpolated_velocity(Utils::Vector3d const &pos) const {
  return std::visit(
      [&](auto &ptr) {
        auto const res = ptr->get_velocity_at_pos(pos * m_conv.pos_to_lb, true);
        assert(res);
        return *res * m_conv.vel_to_md;
      },
      *impl->solver);
}

std::vector<Utils::Vector3d> Solver::get_coupling_interpolated_velocities(
    std::vector<Utils::Vector3d> const &pos) const {
  return std::visit(
      [&](auto &ptr) {
        std::vector<Utils::Vector3d> pos_lb;
        pos_lb.reserve(pos.size());
        for (auto const &pos_md : pos) {
          pos_lb.emplace_back(pos_md * m_conv.pos_to_lb);
        }
        auto res = ptr->get_velocities_at_pos(pos_lb);
        for (auto &v : res) {
          v *= m_conv.vel_to_md;
        }
        return res;
      },
      *impl->solver);
}

void Solver::add_forces_at_pos(std::vector<Utils::Vector3d> const &pos,
                               std::vector<Utils::Vector3d> const &forces) {
  std::visit(
      [&](auto &ptr) {
        std::vector<Utils::Vector3d> pos_lb;
        std::vector<Utils::Vector3d> force_lb;
        pos_lb.reserve(pos.size());
        force_lb.reserve(pos.size());
        for (auto const &pos_md : pos) {
          pos_lb.emplace_back(pos_md * m_conv.pos_to_lb);
        }
        for (auto const &force_md : forces) {
          force_lb.emplace_back(force_md * m_conv.force_to_lb);
        }
        ptr->add_forces_at_pos(pos_lb, force_lb);
      },
      *impl->solver);
}

void Solver::add_force_density(Utils::Vector3d const &pos,
                               Utils::Vector3d const &force_density) {
  std::visit(
      [&](auto &ptr) {
        if (not ptr->add_force_at_pos(pos * m_conv.pos_to_lb,
                                      force_density * m_conv.force_to_lb)) {
          throw std::runtime_error("Cannot apply force to LB");
        }
      },
      *impl->solver);
}

Utils::Vector3d Solver::get_momentum() const {
  check_solver(impl);
  return std::visit([](auto const &ptr) { return ptr->get_momentum(); },
                    *impl->solver);
}

template <> void Solver::set<LBNone>(std::shared_ptr<LBNone> lb_instance) {
  assert(impl);
  assert(not impl->solver.has_value());
  impl->solver = lb_instance;
}

#ifdef WALBERLA
template <>
void Solver::set<LBWalberla>(std::shared_ptr<LBWalberlaBase> lb_fluid,
                             std::shared_ptr<LBWalberlaParams> lb_params) {
  assert(impl);
  assert(not impl->solver.has_value());
  auto const &system = get_system();
  auto lb_instance = std::make_shared<LBWalberla>(lb_fluid, lb_params);
  lb_instance->sanity_checks(system);
  auto const &lebc = system.box_geo->lees_edwards_bc();
  lb_fluid->check_lebc(lebc.shear_direction, lebc.shear_plane_normal);
  impl->solver = lb_instance;
  auto const agrid = lb_instance->get_agrid();
  auto const tau = lb_instance->get_tau();
  m_conv = Conversions{1. / agrid, agrid / tau, tau * tau / agrid};
}
#endif // WALBERLA

} // namespace LB
