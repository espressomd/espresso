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

#pragma once

#include "config/config.hpp"

#ifdef WALBERLA

#include <utils/Vector.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>

class LBWalberlaBase;

namespace LB {

struct LBWalberlaParams {
  LBWalberlaParams(double agrid, double tau) : m_agrid(agrid), m_tau(tau) {}
  double get_agrid() const { return m_agrid; };
  double get_tau() const { return m_tau; };

private:
  double m_agrid;
  double m_tau;
};

struct LBWalberla {
  std::shared_ptr<LBWalberlaBase> lb_fluid;
  std::shared_ptr<LBWalberlaParams> lb_params;
  LBWalberla(std::shared_ptr<LBWalberlaBase> lb_fluid_,
             std::shared_ptr<LBWalberlaParams> lb_params_)
      : lb_fluid{std::move(lb_fluid_)}, lb_params{std::move(lb_params_)} {}
  double get_kT() const;
  auto get_tau() const { return lb_params->get_tau(); }
  auto get_agrid() const { return lb_params->get_agrid(); }
  auto get_lattice_speed() const { return get_agrid() / get_tau(); }
  Utils::VectorXd<9> get_pressure_tensor() const;
  std::optional<Utils::Vector3d>
  get_velocity_at_pos(Utils::Vector3d const &pos,
                      bool consider_points_in_halo) const;
  std::optional<double> get_density_at_pos(Utils::Vector3d const &pos,
                                           bool consider_points_in_halo) const;
  Utils::Vector3d get_momentum() const;
  bool add_force_at_pos(Utils::Vector3d const &pos,
                        Utils::Vector3d const &force);
  void propagate();
  void veto_time_step(double time_step) const;
  void sanity_checks() const;
  void lebc_sanity_checks(unsigned int shear_direction,
                          unsigned int shear_plane_normal) const;

  void on_cell_structure_change() const {}
  void on_boxl_change() const {
    throw std::runtime_error("MD cell geometry change not supported by LB");
  }
  void on_node_grid_change() const {
    throw std::runtime_error("MPI topology change not supported by LB");
  }
  void on_timestep_change() const {
    throw std::runtime_error("Time step change not supported by LB");
  }
  void on_temperature_change() const {
    throw std::runtime_error("Temperature change not supported by LB");
  }
};

} // namespace LB

#endif // WALBERLA
