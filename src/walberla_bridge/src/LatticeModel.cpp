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

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/VTKHandle.hpp>

#include <blockforest/StructuredBlockForest.h>
#include <field/vtk/VTKWriter.h>
#include <vtk/VTKOutput.h>

#include <memory>
#include <sstream>
#include <string>

std::shared_ptr<VTKHandle> LatticeModel::create_vtk(
    int delta_N, int initial_count, int flag_observables,
    units_map const &units_conversion, std::string const &identifier,
    std::string const &base_folder, std::string const &prefix) {

  using walberla::uint_c;

  // VTKOutput object must be unique
  std::stringstream unique_identifier;
  unique_identifier << base_folder << "/" << identifier;
  std::string const vtk_uid = unique_identifier.str();
  if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end() or
      m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
    throw vtk_runtime_error(vtk_uid, "already exists");
  }

  // instantiate VTKOutput object
  auto const &blocks = get_lattice().get_blocks();
  auto const write_freq = (delta_N) ? static_cast<unsigned int>(delta_N) : 1u;
  auto vtk_obj = walberla::vtk::createVTKOutput_BlockData(
      blocks, identifier, uint_c(write_freq), uint_c(0), false, base_folder,
      prefix, true, true, true, true, uint_c(initial_count));

  // add filters
  register_vtk_field_filters(*vtk_obj);

  // add writers
  register_vtk_field_writers(*vtk_obj, units_conversion, flag_observables);

  auto vtk_handle = std::make_shared<VTKHandle>(vtk_obj, initial_count, true);
  if (delta_N) {
    m_vtk_auto[vtk_uid] = vtk_handle;
  } else {
    m_vtk_manual[vtk_uid] = vtk_handle;
  }
  return vtk_handle;
}

void LatticeModel::write_vtk(std::string const &vtk_uid) {
  if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end()) {
    throw vtk_runtime_error(vtk_uid, "is an automatic observable");
  }
  if (m_vtk_manual.find(vtk_uid) == m_vtk_manual.end()) {
    throw vtk_runtime_error(vtk_uid, "doesn't exist");
  }
  auto &vtk_handle = m_vtk_manual[vtk_uid];
  walberla::vtk::writeFiles(vtk_handle->ptr)();
  vtk_handle->execution_count++;
}

void LatticeModel::switch_vtk(std::string const &vtk_uid, bool status) {
  if (m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
    throw vtk_runtime_error(vtk_uid, "is a manual observable");
  }
  if (m_vtk_auto.find(vtk_uid) == m_vtk_auto.end()) {
    throw vtk_runtime_error(vtk_uid, "doesn't exist");
  }
  m_vtk_auto[vtk_uid]->enabled = status;
}
