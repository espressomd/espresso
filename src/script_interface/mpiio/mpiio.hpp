/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "script_interface/ObjectHandle.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/system/System.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/io/mpiio/mpiio.hpp"
#include "core/system/System.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace MPIIO {

class MPIIOScript : public ObjectHandle {
  std::weak_ptr<System::System> m_system;
  std::shared_ptr<Mpiio::write_buffers> m_buffers;

public:
  void do_construct(VariantMap const &parameters) override {
    m_system = get_value<std::shared_ptr<System::System>>(parameters, "system");
    m_buffers = std::make_shared<Mpiio::write_buffers>();
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &parameters) override {

    auto prefix = get_value<std::string>(parameters.at("prefix"));
    auto pos = get_value<bool>(parameters.at("pos"));
    auto vel = get_value<bool>(parameters.at("vel"));
    auto typ = get_value<bool>(parameters.at("typ"));
    auto bnd = get_value<bool>(parameters.at("bond"));

    auto const fields = ((pos) ? Mpiio::MPIIO_OUT_POS : Mpiio::MPIIO_OUT_NON) |
                        ((vel) ? Mpiio::MPIIO_OUT_VEL : Mpiio::MPIIO_OUT_NON) |
                        ((typ) ? Mpiio::MPIIO_OUT_TYP : Mpiio::MPIIO_OUT_NON) |
                        ((bnd) ? Mpiio::MPIIO_OUT_BND : Mpiio::MPIIO_OUT_NON);

    if (name == "write") {
      auto const system_si = m_system.lock();
      auto &system = system_si->get_system();
      auto &cell_structure = *system.cell_structure;
      auto &bonded_ias = *system.bonded_ias;
      Mpiio::mpi_mpiio_common_write(prefix, fields, bonded_ias,
                                    cell_structure.local_particles(),
                                    *m_buffers);
    } else if (name == "read") {
      auto const system_si = m_system.lock();
      auto &system = system_si->get_system();
      auto &cell_structure = *system.cell_structure;
      Mpiio::mpi_mpiio_common_read(prefix, fields, cell_structure);
    }

    return {};
  }
};

} // namespace MPIIO
} // namespace ScriptInterface
