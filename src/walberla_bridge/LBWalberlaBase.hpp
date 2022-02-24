/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#ifndef LB_WALBERLA_BASE_HPP
#define LB_WALBERLA_BASE_HPP

/**
 * @file
 * @ref LBWalberlaBase provides the public interface of the LB
 * waLBerla bridge. It relies on type erasure to hide the waLBerla
 * implementation details from the ESPResSo core. It is implemented
 * by @ref walberla::LBWalberlaImpl.
 */

#include "LatticeWalberla.hpp"
#include "LeesEdwardsPack.hpp"
#include "VTKHandle.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

/** Class that runs and controls the LB on WaLBerla.  */
class LBWalberlaBase {
public:
  virtual ~LBWalberlaBase() = default;

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

  /** @brief Configure an unthermalized collision model. */
  virtual void set_collision_model() = 0;

  /** @brief Configure a thermalized collision model. */
  virtual void set_collision_model(double kT, unsigned int seed) = 0;

  /** @brief Configure a thermalized collision model for Lees-Edwards. */
  virtual void
  set_collision_model(std::unique_ptr<LeesEdwardsPack> &&lees_edwards_pack) = 0;

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

  /** @brief Get the underlying lattice. */
  virtual LatticeWalberla const &lattice() const = 0;

  /** @brief Create a VTK observable.
   *
   *  @param delta_N          Write frequency, if 0 write a single frame,
   *                          otherwise add a callback to write every
   *                          @p delta_N LB steps to a new file
   *  @param initial_count    Initial execution count
   *  @param flag_observables Which observables to measure (OR'ing of
   *                          @ref OutputVTK values)
   *  @param identifier       Name of the VTK dataset
   *  @param base_folder      Path to the VTK folder
   *  @param prefix           Prefix of the VTK files
   */
  virtual std::shared_ptr<VTKHandle> create_vtk(int delta_N, int initial_count,
                                                int flag_observables,
                                                std::string const &identifier,
                                                std::string const &base_folder,
                                                std::string const &prefix) = 0;
  /** @brief Write a VTK observable to disk.
   *
   *  @param vtk_uid          Name of the VTK object
   */
  virtual void write_vtk(std::string const &vtk_uid) = 0;

  /** @brief Toggle a VTK observable on/off.
   *
   *  @param vtk_uid          Name of the VTK object
   *  @param status           @c true to switch on, @c false to switch off
   */
  virtual void switch_vtk(std::string const &vtk_uid, bool status) = 0;
};

#endif // LB_WALBERLA_BASE_HPP
