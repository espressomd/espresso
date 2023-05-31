/*
 * Copyright (C) 2020-2023 The ESPResSo project
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

#include <walberla_bridge/utils/ResourceManager.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace walberla::vtk {
// forward declare
class VTKOutput;
} // namespace walberla::vtk

/** @brief Handle to a VTK object */
class VTKHandle {
  std::unique_ptr<ResourceManager> m_vtk_resources_lock;

public:
  VTKHandle(std::shared_ptr<walberla::vtk::VTKOutput> sp, int ec, bool en)
      : ptr(std::move(sp)), execution_count(ec), enabled(en) {
    m_vtk_resources_lock = ::walberla::get_vtk_dependent_resources();
  }
  ~VTKHandle() {
    // vtk objects must be cleared *before* the MPI resources can be freed,
    // because file handles need to be closed on all ranks
    ptr.reset();
    m_vtk_resources_lock.reset();
  }

  std::shared_ptr<walberla::vtk::VTKOutput> ptr;
  int execution_count;
  bool enabled;
};

/** @brief LB statistics to write to VTK files */
enum class OutputVTK : int {
  density = 1 << 0,
  velocity_vector = 1 << 1,
  pressure_tensor = 1 << 2,
};

/** @brief EK statistics to write to VTK files */
enum class EKOutputVTK : int {
  density = 1 << 0,
};

class vtk_runtime_error : public std::runtime_error {
public:
  explicit vtk_runtime_error(std::string const &vtk_uid,
                             std::string const &reason)
      : std::runtime_error("VTKOutput object '" + vtk_uid + "' " + reason) {}
};
