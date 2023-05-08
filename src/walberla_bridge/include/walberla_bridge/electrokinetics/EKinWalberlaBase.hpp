/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include <walberla_bridge/LatticeModel.hpp>

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cstddef>
#include <vector>

/** @brief Interface of a lattice-based electrokinetic model. */
class EKinWalberlaBase : public LatticeModel {
public:
  /** @brief Integrate EKin for one time step */
  virtual void integrate(std::size_t potential_id, std::size_t velocity_id,
                         std::size_t force_id) = 0;

  /** @brief perform ghost communication of densities */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized fluxes */
  [[nodiscard]] virtual std::size_t stencil_size() const noexcept = 0;

  /** @brief Set node density. */
  virtual bool set_node_density(Utils::Vector3i const &node,
                                double density) = 0;

  /** @brief Get node density. */
  [[nodiscard]] virtual boost::optional<double>
  get_node_density(Utils::Vector3i const &node,
                   bool consider_ghosts = false) const = 0;

  /** @brief Set slice density. */
  virtual void set_slice_density(Utils::Vector3i const &lower_corner,
                                 Utils::Vector3i const &upper_corner,
                                 std::vector<double> const &density) = 0;

  /** @brief Get slice density. */
  [[nodiscard]] virtual std::vector<double>
  get_slice_density(Utils::Vector3i const &lower_corner,
                    Utils::Vector3i const &upper_corner) const = 0;

  /** @brief Set node flux boundary conditions. */
  virtual bool set_node_flux_boundary(Utils::Vector3i const &node,
                                      Utils::Vector3d const &flux) = 0;

  /** @brief Get node flux boundary conditions. */
  [[nodiscard]] virtual boost::optional<Utils::Vector3d>
  get_node_flux_at_boundary(Utils::Vector3i const &node,
                            bool consider_ghosts = false) const = 0;

  /** @brief Set slice flux boundary conditions. */
  virtual void set_slice_flux_boundary(
      Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
      std::vector<boost::optional<Utils::Vector3d>> const &flux) = 0;

  /** @brief Get slice flux boundary conditions. */
  [[nodiscard]] virtual std::vector<boost::optional<Utils::Vector3d>>
  get_slice_flux_at_boundary(Utils::Vector3i const &lower_corner,
                             Utils::Vector3i const &upper_corner) const = 0;

  virtual bool remove_node_from_flux_boundary(Utils::Vector3i const &node) = 0;

  /** @brief Set node density boundary conditions. */
  virtual bool set_node_density_boundary(Utils::Vector3i const &node,
                                         double density) = 0;

  /** @brief Get node density boundary conditions. */
  [[nodiscard]] virtual boost::optional<double>
  get_node_density_at_boundary(Utils::Vector3i const &node,
                               bool consider_ghosts = false) const = 0;

  /** @brief Set slice density boundary conditions. */
  virtual void set_slice_density_boundary(
      Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
      std::vector<boost::optional<double>> const &density) = 0;

  /** @brief Get slice density boundary conditions. */
  [[nodiscard]] virtual std::vector<boost::optional<double>>
  get_slice_density_at_boundary(Utils::Vector3i const &lower_corner,
                                Utils::Vector3i const &upper_corner) const = 0;

  virtual bool
  remove_node_from_density_boundary(Utils::Vector3i const &node) = 0;

  /** @brief Check if node has flux boundary conditions. */
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_flux_boundary(Utils::Vector3i const &node,
                            bool consider_ghosts = false) const = 0;

  /** @brief Check if node has density boundary conditions. */
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_density_boundary(Utils::Vector3i const &node,
                               bool consider_ghosts = false) const = 0;

  /** @brief Check if node has any boundary conditions. */
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &node,
                       bool consider_ghosts = false) const = 0;

  /** @brief Check if slice has any boundary conditions. */
  [[nodiscard]] virtual std::vector<bool>
  get_slice_is_boundary(Utils::Vector3i const &lower_corner,
                        Utils::Vector3i const &upper_corner) const = 0;

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

  virtual void set_diffusion(double diffusion) = 0;
  virtual void set_kT(double kT) = 0;
  virtual void set_valency(double valency) = 0;
  virtual void set_advection(bool advection) = 0;
  virtual void set_friction_coupling(bool friction_coupling) = 0;
  virtual void set_ext_efield(Utils::Vector3d const &field) = 0;

  [[nodiscard]] virtual std::size_t get_density_id() const noexcept = 0;

  ~EKinWalberlaBase() override = default;
};
