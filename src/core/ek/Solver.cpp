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

#include "ek/Implementation.hpp"
#include "ek/Solver.hpp"
#include "ek/utils.hpp"

#include "ek/EKNone.hpp"
#include "ek/EKWalberla.hpp"

#include "integrate.hpp"
#include "system/System.hpp"

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

namespace EK {

Solver::Solver() { impl = std::make_unique<Implementation>(); }

static auto is_solver_set(std::unique_ptr<Solver::Implementation> const &ptr) {
  return ptr != nullptr and ptr->solver.has_value();
}

static void check_solver(std::unique_ptr<Solver::Implementation> const &ptr) {
  if (not is_solver_set(ptr)) {
    throw NoEKActive{};
  }
}

bool Solver::is_solver_set() const { return EK::is_solver_set(impl); }

void Solver::reset() { System::get_system().ek.impl->solver = std::nullopt; }

bool Solver::is_ready_for_propagation() const {
  return is_solver_set() and
         std::visit([](auto &ptr) { return ptr->is_ready_for_propagation(); },
                    *impl->solver);
}

void Solver::propagate() {
  check_solver(impl);
  std::visit([](auto &ptr) { ptr->propagate(); }, *impl->solver);
}

void Solver::sanity_checks() const {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->sanity_checks(); }, *impl->solver);
  }
}

void Solver::veto_time_step(double time_step) const {
  if (impl->solver) {
    std::visit([=](auto &ptr) { ptr->veto_time_step(time_step); },
               *impl->solver);
  }
}

void Solver::on_cell_structure_change() {
  if (impl->solver) {
    auto &solver = *impl->solver;
    std::visit([](auto &ptr) { ptr->on_cell_structure_change(); }, solver);
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

double Solver::get_tau() const {
  check_solver(impl);
  return std::visit([](auto &ptr) { return ptr->get_tau(); }, *impl->solver);
}

template <> void Solver::set<EKNone>(std::shared_ptr<EKNone> ek_instance) {
  assert(impl);
  assert(not impl->solver.has_value());
  impl->solver = ek_instance;
}

#ifdef WALBERLA
template <>
void Solver::set<EKWalberla>(std::shared_ptr<EKWalberla> ek_instance) {
  assert(impl);
  assert(not impl->solver.has_value());
  ek_instance->sanity_checks();
  impl->solver = ek_instance;
}
#endif // WALBERLA

} // namespace EK
