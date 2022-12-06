/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

/**
 * @file
 * @ref LBWalberlaBase provides the public interface of the LB
 * waLBerla bridge. It relies on type erasure to hide the waLBerla
 * implementation details from the ESPResSo core. It is implemented
 * by @ref walberla::LBWalberlaImpl.
 */

#include "LeesEdwardsPack.hpp"
#include "walberla_bridge/LatticeModel.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

/** @brief Interface of a lattice-based fluid model. */
class LBWalberlaBase : public LatticeModel {
public:
  ~LBWalberlaBase() override = default;

  /** @brief Integrate LB for one time step. */
  virtual void integrate() = 0;

  /** @brief Perform ghost communication of PDF and applied forces. */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized velocities in the PDF. */
  virtual std::size_t stencil_size() const = 0;

  /** @brief Get node velocity. */
  virtual boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i &node,
                    bool consider_ghosts = false) const = 0;

  /** @brief Set node velocity. */
  virtual bool set_node_velocity(const Utils::Vector3i &node,
                                 const Utils::Vector3d &v) = 0;

  /** @brief Get interpolated velocities at a position. */
  virtual boost::optional<Utils::Vector3d>
  get_velocity_at_pos(const Utils::Vector3d &position,
                      bool consider_points_in_halo = false) const = 0;

  /** @brief Get interpolated densities at a position. */
  virtual boost::optional<double> get_interpolated_density_at_pos(
      const Utils::Vector3d &position,
      bool consider_points_in_halo = false) const = 0;

  /**
   * @brief Interpolate a force to the stored forces to be applied on nodes
   * in the next time step.
   */
  virtual bool add_force_at_pos(const Utils::Vector3d &position,
                                const Utils::Vector3d &force) = 0;

  /** @brief Get stored force to be applied on node in the next time step. */
  virtual boost::optional<Utils::Vector3d>
  get_node_force_to_be_applied(const Utils::Vector3i &node) const = 0;

  /** @brief Get stored force that was applied on node in the last time step. */
  virtual boost::optional<Utils::Vector3d>
  get_node_last_applied_force(const Utils::Vector3i &node,
                              bool consider_ghosts = false) const = 0;

  /** @brief Set stored force that was applied on node in the last time step. */
  virtual bool set_node_last_applied_force(const Utils::Vector3i &node,
                                           const Utils::Vector3d &force) = 0;

  /** @brief Get node populations. */
  virtual boost::optional<std::vector<double>>
  get_node_pop(const Utils::Vector3i &node,
               bool consider_ghosts = false) const = 0;

  /** @brief Set node populations. */
  virtual bool set_node_pop(const Utils::Vector3i &node,
                            std::vector<double> const &population) = 0;

  /** @brief Set node density. */
  virtual bool set_node_density(const Utils::Vector3i &node,
                                double density) = 0;

  /** @brief Get node density. */
  virtual boost::optional<double>
  get_node_density(const Utils::Vector3i &node,
                   bool consider_ghosts = false) const = 0;

  /** @brief Get node velocity boundary conditions. */
  virtual boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const = 0;

  /** @brief Set node velocity boundary conditions. */
  virtual bool set_node_velocity_at_boundary(const Utils::Vector3i &node,
                                             const Utils::Vector3d &v) = 0;

  /** @brief Get (stored) force applied on node due to boundary condition. */
  virtual boost::optional<Utils::Vector3d>
  get_node_boundary_force(const Utils::Vector3i &node) const = 0;

  /** @brief Remove a node from the boundaries. */
  virtual bool remove_node_from_boundary(const Utils::Vector3i &node) = 0;

  /** @brief Check if node has velocity boundary conditions. */
  virtual boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const = 0;

  /** @brief Rebuild the UBB field. This is an expensive operation. */
  virtual void reallocate_ubb_field() = 0;

  /** @brief Clear the boundary flag field and the UBB field. */
  virtual void clear_boundaries() = 0;

  /** @brief Update boundary conditions from a rasterized shape. */
  virtual void update_boundary_from_shape(std::vector<int> const &,
                                          std::vector<double> const &) = 0;

  /** @brief Configure the default collision model. */
  virtual void set_collision_model(double kT, unsigned int seed) = 0;

  /** @brief Configure a thermalized collision model for Lees-Edwards. */
  virtual void
  set_collision_model(std::unique_ptr<LeesEdwardsPack> &&lees_edwards_pack) = 0;

  /** @brief Check Lees-Edwards boundary conditions. */
  virtual void check_lebc(int shear_direction,
                          int shear_plane_normal) const = 0;

  /** @brief Get node pressure tensor. */
  virtual boost::optional<Utils::VectorXd<9>>
  get_node_pressure_tensor(const Utils::Vector3i &node) const = 0;

  /** @brief Calculate average pressure tensor of the local domain. */
  virtual Utils::VectorXd<9> get_pressure_tensor() const = 0;

  /** @brief Calculate momentum of the local domain. */
  virtual Utils::Vector3d get_momentum() const = 0;

  /** @brief Set a global external force. */
  virtual void set_external_force(const Utils::Vector3d &ext_force) = 0;

  /** @brief Get the global external force. */
  virtual Utils::Vector3d get_external_force() const = 0;

  /** @brief Set the fluid viscosity. */
  virtual void set_viscosity(double viscosity) = 0;

  /** @brief Get the fluid viscosity. */
  virtual double get_viscosity() const = 0;

  /** @brief Get the fluid density. */
  virtual double get_density() const = 0;

  /** @brief Get the fluid temperature (if thermalized). */
  virtual double get_kT() const = 0;

  /** @brief Set the RNG counter (if thermalized). */
  virtual uint64_t get_rng_state() const = 0;

  /** @brief Set the rng state of thermalized LBs */
  virtual void set_rng_state(uint64_t counter) = 0;

  /** @brief get the velocity field id */
  [[nodiscard]] virtual std::size_t get_velocity_field_id() const = 0;

  /** @brief get the force field id */
  [[nodiscard]] virtual std::size_t get_force_field_id() const = 0;
};
