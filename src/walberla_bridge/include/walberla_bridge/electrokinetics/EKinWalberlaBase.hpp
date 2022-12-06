/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "walberla_bridge/LatticeModel.hpp"

#include <boost/optional.hpp>

#include <utils/Vector.hpp>

#include <cstddef>
#include <string>

/** @brief Interface of a lattice-based electrokinetic model. */
class EKinWalberlaBase : public LatticeModel {
public:
  /** @brief Integrate EKin for one time step */
  virtual void integrate(const std::size_t &potential_id,
                         const std::size_t &velocity_id,
                         const std::size_t &force_id) = 0;

  /** @brief perform ghost communication of densities */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized fluxes */
  [[nodiscard]] virtual std::size_t stencil_size() const = 0;

  // Density
  virtual bool set_node_density(const Utils::Vector3i &node,
                                double density) = 0;
  [[nodiscard]] virtual boost::optional<double>
  get_node_density(const Utils::Vector3i &node) const = 0;

  virtual bool set_node_flux_boundary(const Utils::Vector3i &node,
                                      const Utils::Vector3d &flux) = 0;
  [[nodiscard]] virtual boost::optional<Utils::Vector3d>
  get_node_flux_at_boundary(const Utils::Vector3i &node) const = 0;
  virtual bool remove_node_from_flux_boundary(const Utils::Vector3i &node) = 0;
  virtual bool set_node_density_boundary(const Utils::Vector3i &node,
                                         double density) = 0;
  [[nodiscard]] virtual boost::optional<double>
  get_node_density_at_boundary(const Utils::Vector3i &node) const = 0;
  virtual bool
  remove_node_from_density_boundary(const Utils::Vector3i &node) = 0;
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_flux_boundary(const Utils::Vector3i &node,
                            bool consider_ghosts = false) const = 0;
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_density_boundary(const Utils::Vector3i &node,
                               bool consider_ghosts = false) const = 0;
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const = 0;
  virtual void clear_flux_boundaries() = 0;
  virtual void clear_density_boundaries() = 0;

  virtual void update_flux_boundary_from_shape(std::vector<int> const &,
                                               std::vector<double> const &) = 0;
  virtual void
  update_density_boundary_from_shape(std::vector<int> const &,
                                     std::vector<double> const &) = 0;

  // Global parameters
  [[nodiscard]] virtual double get_diffusion() const noexcept = 0;
  [[nodiscard]] virtual double get_kT() const noexcept = 0;
  [[nodiscard]] virtual double get_valency() const noexcept = 0;
  [[nodiscard]] virtual bool get_advection() const noexcept = 0;
  [[nodiscard]] virtual bool get_friction_coupling() const noexcept = 0;
  [[nodiscard]] virtual Utils::Vector3d get_ext_efield() const noexcept = 0;
  [[nodiscard]] virtual bool is_double_precision() const noexcept = 0;

  virtual void set_diffusion(double diffusion) noexcept = 0;
  virtual void set_kT(double kT) noexcept = 0;
  virtual void set_valency(double valency) noexcept = 0;
  virtual void set_advection(bool advection) noexcept = 0;
  virtual void set_friction_coupling(bool friction_coupling) noexcept = 0;
  virtual void set_ext_efield(const Utils::Vector3d &field) noexcept = 0;

  [[nodiscard]] virtual std::size_t get_density_id() const noexcept = 0;

  ~EKinWalberlaBase() override = default;
};
