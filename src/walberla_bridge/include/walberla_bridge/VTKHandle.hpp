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

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace walberla::vtk {
// forward declare
class VTKOutput;
} // namespace walberla::vtk

/** @brief Handle to a VTK object */
struct VTKHandle {
  VTKHandle(std::shared_ptr<walberla::vtk::VTKOutput> sp, int ec, bool en)
      : ptr(std::move(sp)), execution_count(ec), enabled(en) {}
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
