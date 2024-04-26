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

#include "system/Leaf.hpp"

#include "utils.hpp"

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <memory>
#include <optional>

namespace LB {

struct Solver : public System::Leaf<Solver> {
  struct Implementation;

  Solver();

  /** @brief Return true if a LB solver is active. */
  [[nodiscard]] bool is_solver_set() const;

  /** @brief Remove the LB solver. */
  void reset();

  /**
   * @brief Set the LB solver.
   * For developers: a specialization must exist for every new LB type.
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
   * @brief Propagate the LB fluid.
   */
  void propagate();

  /**
   * @brief Perform a full initialization of the lattice-Boltzmann system.
   * All derived parameters and the fluid are reset to their default values.
   */
  void init() const {}

  /**
   * @brief Perform LB parameter and boundary velocity checks.
   */
  void sanity_checks() const;

  /**
   * @brief Check if a MD time step is compatible with the LB tau.
   */
  void veto_time_step(double time_step) const;

  /**
   * @brief Check if a thermostat is compatible with the LB temperature.
   */
  void veto_kT(double kT) const;

  /**
   * @brief Perform LB LEbc parameter checks.
   */
  void lebc_sanity_checks(unsigned int shear_direction,
                          unsigned int shear_plane_normal) const;

  /**
   * @brief Get the LB time step.
   */
  double get_tau() const;

  /**
   * @brief Get the LB grid spacing.
   */
  double get_agrid() const;

  /**
   * @brief Get the thermal energy parameter of the LB fluid.
   */
  double get_kT() const;

  /**
   * @brief Get the lattice speed (agrid/tau).
   */
  auto get_lattice_speed() const { return get_agrid() / get_tau(); }

  Utils::VectorXd<9> get_pressure_tensor() const;

  Utils::Vector3d get_momentum() const;

  /**
   * @brief Calculate the interpolated fluid velocity in LB units.
   * Use this function in MPI-parallel code. The LB ghost layer is ignored.
   * @param pos Position in MD units at which the velocity is to be calculated.
   * @retval interpolated fluid velocity.
   */
  std::optional<Utils::Vector3d>
  get_interpolated_velocity(Utils::Vector3d const &pos) const;

  /**
   * @brief Calculate the interpolated fluid density in LB units.
   * Use this function in MPI-parallel code. The LB ghost layer is ignored.
   * @param pos Position in MD units at which the density is to be calculated.
   * @retval interpolated fluid density.
   */
  std::optional<double>
  get_interpolated_density(Utils::Vector3d const &pos) const;

  /**
   * @brief Calculate the interpolated fluid velocity in MD units.
   * Special method used only for particle coupling. Uses the LB ghost layer.
   * @param pos Position in MD units at which the velocity is to be calculated.
   * @retval interpolated fluid velocity.
   */
  Utils::Vector3d
  get_coupling_interpolated_velocity(Utils::Vector3d const &pos) const;

  /**
   * @brief Add a force density to the fluid at the given position.
   * @param pos            Position at which the force density is to be applied.
   * @param force_density  Force density to apply.
   */
  void add_force_density(Utils::Vector3d const &pos,
                         Utils::Vector3d const &force_density);

  void on_boxl_change();
  void on_node_grid_change();
  void on_cell_structure_change();
  void on_timestep_change();
  void on_temperature_change();
  void veto_boxl_change() const;

private:
  /** @brief Pointer-to-implementation. */
  std::unique_ptr<Implementation> impl;
};

} // namespace LB
