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

#include "../../VTKHandle.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"

#include <map>
#include <memory>
#include <string>

/** @brief Abstract representation of a lattice-based model. */
class LatticeModel {
protected:
  /** VTK writers that are executed automatically */
  std::map<std::string, std::shared_ptr<VTKHandle>> m_vtk_auto;
  /** VTK writers that are executed manually */
  std::map<std::string, std::shared_ptr<VTKHandle>> m_vtk_manual;

  virtual void register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                                          int flag_observables) = 0;

  virtual void
  register_vtk_field_filters(walberla::vtk::VTKOutput &vtk_obj) = 0;

  virtual void integrate_vtk_writers() = 0;

public:
  virtual ~LatticeModel() = default;

  /** @brief Get the underlying lattice. */
  virtual LatticeWalberla const &get_lattice() const noexcept = 0;

  /** @brief Create a VTK observable.
   *
   *  @param delta_N          Write frequency, if 0 write a single frame,
   *                          otherwise add a callback to write every
   *                          @p delta_N EK steps to a new file
   *  @param initial_count    Initial execution count
   *  @param flag_observables Which observables to measure (OR'ing of
   *                          @ref OutputVTK  or @ref EKOutputVTK values)
   *  @param identifier       Name of the VTK dataset
   *  @param base_folder      Path to the VTK folder
   *  @param prefix           Prefix of the VTK files
   */
  std::shared_ptr<VTKHandle> create_vtk(int delta_N, int initial_count,
                                        int flag_observables,
                                        std::string const &identifier,
                                        std::string const &base_folder,
                                        std::string const &prefix);

  /** @brief Write a VTK observable to disk.
   *
   *  @param vtk_uid          Name of the VTK object
   */
  void write_vtk(std::string const &vtk_uid);

  /** @brief Toggle a VTK observable on/off.
   *
   *  @param vtk_uid          Name of the VTK object
   *  @param status           @c true to switch on, @c false to switch off
   */
  void switch_vtk(std::string const &vtk_uid, bool status);
};
