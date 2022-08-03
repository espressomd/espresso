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

#ifndef ESPRESSO_EKINWALBERLABASE_HPP
#define ESPRESSO_EKINWALBERLABASE_HPP

#include "LatticeWalberla.hpp"

#include <boost/optional.hpp>
#include <string>

#include "utils/Vector.hpp"

#include "VTKHandle.hpp"

#include <domain_decomposition/BlockDataID.h>

/** Class that runs and controls the EK on WaLBerla
 */
class EKinWalberlaBase {
public:
  /** @brief Integrate EKin for one time step */
  virtual void integrate(const std::size_t &potential_id,
                         const std::size_t &velocity_id,
                         const std::size_t &force_id) = 0;

  /** @brief perform ghost communication of densities */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized fluxes */
  [[nodiscard]] virtual size_t stencil_size() const = 0;

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

  [[nodiscard]] virtual walberla::BlockDataID get_density_id() const
      noexcept = 0;

  [[nodiscard]] virtual LatticeWalberla &get_lattice() const noexcept = 0;

  /** @brief Create a VTK observable.
   *
   *  @param delta_N          Write frequency, if 0 write a single frame,
   *                          otherwise add a callback to write every
   *                          @p delta_N EK steps to a new file
   *  @param initial_count    Initial execution count
   *  @param flag_observables Which observables to measure (OR'ing of
   *                          @ref EKOutputVTK values)
   *  @param identifier       Name of the VTK dataset
   *  @param base_folder      Path to the VTK folder
   *  @param prefix           Prefix of the VTK files
   */
  [[nodiscard]] virtual std::shared_ptr<VTKHandle>
  create_vtk(int delta_N, int initial_count, int flag_observables,
             std::string const &identifier, std::string const &base_folder,
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

  virtual ~EKinWalberlaBase() = default;
};

#endif // ESPRESSO_EKINWALBERLABASE_HPP
