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

#include "utils.hpp"

namespace System {
class System;
}

namespace LB {

struct LBNone {
  void propagate() { throw NoLBActive{}; }
  double get_agrid() const { throw NoLBActive{}; }
  double get_tau() const { throw NoLBActive{}; }
  double get_kT() const { throw NoLBActive{}; }
  Utils::VectorXd<9> get_pressure_tensor() const { throw NoLBActive{}; }
  std::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d const &,
                                                     bool) const {
    throw NoLBActive{};
  }
  std::optional<double> get_density_at_pos(Utils::Vector3d const &,
                                           bool) const {
    throw NoLBActive{};
  }
  bool add_force_at_pos(Utils::Vector3d const &pos,
                        Utils::Vector3d const &force) const {
    throw NoLBActive{};
  }
  Utils::Vector3d get_momentum() const { throw NoLBActive{}; }
  void veto_time_step(double) const { throw NoLBActive{}; }
  void veto_kT(double) const { throw NoLBActive{}; }
  void sanity_checks(System::System const &) const { throw NoLBActive{}; }
  void lebc_sanity_checks(unsigned int, unsigned int) const {
    throw NoLBActive{};
  }
  void veto_boxl_change() const { throw NoLBActive{}; }
  void on_cell_structure_change() const { throw NoLBActive{}; }
  void on_boxl_change() const { throw NoLBActive{}; }
  void on_node_grid_change() const { throw NoLBActive{}; }
  void on_timestep_change() const { throw NoLBActive{}; }
  void on_temperature_change() const { throw NoLBActive{}; }
};

} // namespace LB
