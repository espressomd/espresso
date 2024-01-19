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

namespace EK {

struct EKNone {
  bool is_ready_for_propagation() const { throw NoEKActive{}; }
  void propagate() { throw NoEKActive{}; }
  double get_tau() const { throw NoEKActive{}; }
  void veto_time_step(double) const { throw NoEKActive{}; }
  void veto_kT(double) const { throw NoEKActive{}; }
  void sanity_checks(System::System const &) const { throw NoEKActive{}; }
  void on_cell_structure_change() const { throw NoEKActive{}; }
  void on_boxl_change() const { throw NoEKActive{}; }
  void on_node_grid_change() const { throw NoEKActive{}; }
  void on_timestep_change() const { throw NoEKActive{}; }
  void on_temperature_change() const { throw NoEKActive{}; }
};

} // namespace EK
