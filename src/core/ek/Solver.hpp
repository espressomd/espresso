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

#include <cassert>
#include <cmath>
#include <memory>
#include <optional>

namespace EK {

struct Solver {
  struct Implementation;

  Solver();

  /** @brief Return true if an @c EK solver is active. */
  [[nodiscard]] bool is_solver_set() const;

  /** @brief Return true if an @c EK solver can be propagated. */
  bool is_ready_for_propagation() const;

  /** @brief Remove the @c EK solver. */
  void reset();

  /**
   * @brief Set the @c EK solver.
   * For developers: a specialization must exist for every new @c EK type.
   */
  template <typename LB, class... Args> void set(Args... args);

  /**
   * @brief Connector to the implementation internal state.
   * For developers: use this mechanism to access the underlying variant.
   */
  template <class Connector> void connect(Connector &&connector) const {
    assert(impl != nullptr);
    connector(*this->impl);
  }

  /**
   * @brief Propagate the @c EK species.
   */
  void propagate();

  /**
   * @brief Perform a full initialization of the lattice-Boltzmann system.
   * All derived parameters and the fluid are reset to their default values.
   */
  void init() const {}

  /**
   * @brief Get the @c EK time step.
   */
  double get_tau() const;

  auto get_steps_per_md_step(double time_step) const {
    return static_cast<int>(std::round(get_tau() / time_step));
  }

  /**
   * @brief Perform @c EK parameter checks.
   */
  void sanity_checks() const;

  /**
   * @brief Check if a MD time step is compatible with the @c EK tau.
   */
  void veto_time_step(double time_step) const;

  void on_boxl_change();
  void on_node_grid_change();
  void on_cell_structure_change();
  void on_timestep_change();
  void on_temperature_change();

private:
  /** @brief Pointer-to-implementation. */
  std::unique_ptr<Implementation> impl;
};

} // namespace EK
