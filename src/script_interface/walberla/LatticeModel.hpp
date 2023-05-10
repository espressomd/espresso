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

#include "LatticeWalberla.hpp"

#include <walberla_bridge/LatticeModel.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

template <class Method, class VTKHandle>
class LatticeModel : public AutoParameters<LatticeModel<Method, VTKHandle>> {
protected:
  std::shared_ptr<LatticeWalberla> m_lattice;
  std::shared_ptr<Method> m_instance;
  std::vector<std::shared_ptr<VTKHandle>> m_vtk_writers;

  virtual ::LatticeModel::units_map
  get_latice_to_md_units_conversion() const = 0;

  auto find_vtk(std::shared_ptr<VTKHandle> const &vtk) const {
    return std::find(m_vtk_writers.begin(), m_vtk_writers.end(), vtk);
  }

  auto serialize_vtk_writers() const {
    return make_vector_of_variants(m_vtk_writers);
  }

public:
  Variant do_call_method(std::string const &method_name,
                         VariantMap const &params) override {
    if (method_name == "add_vtk_writer") {
      auto vtk = get_value<std::shared_ptr<VTKHandle>>(params, "vtk");
      auto const needle = find_vtk(vtk);
      ObjectHandle::context()->parallel_try_catch([&]() {
        if (needle != m_vtk_writers.end()) {
          throw std::runtime_error(
              "VTK object is already attached to this lattice");
        }
        vtk->attach_to_lattice(m_instance, get_latice_to_md_units_conversion());
        m_vtk_writers.emplace_back(vtk);
      });
      return {};
    }
    if (method_name == "remove_vtk_writer") {
      auto const vtk = get_value<std::shared_ptr<VTKHandle>>(params, "vtk");
      auto const needle = find_vtk(vtk);
      ObjectHandle::context()->parallel_try_catch([&]() {
        if (needle == m_vtk_writers.end()) {
          throw std::runtime_error(
              "VTK object is not attached to this lattice");
        }
        vtk->detach_from_lattice();
      });
      m_vtk_writers.erase(needle);
      return {};
    }
    if (method_name == "clear_vtk_writers") {
      for (auto const &vtk : m_vtk_writers) {
        vtk->detach_from_lattice();
      }
      m_vtk_writers.clear();
    }
    return {};
  }
};

} // namespace ScriptInterface::walberla
